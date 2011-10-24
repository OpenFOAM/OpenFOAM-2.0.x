/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "codedBaseClass.H"
#include "dictionary.H"
#include "SHA1Digest.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"
#include "stringOps.H"
#include "dlLibraryTable.H"
#include "Time.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void* Foam::codedBaseClass::loadLibrary
(
    const fileName& libPath,
    const string& globalFuncName,
    const dictionary& contextDict
) const
{
    void* lib = 0;

    // avoid compilation by loading an existing library
    if (!libPath.empty())
    {
        dlLibraryTable& libs = this->libs();

        if (libs.open(libPath, false))
        {
            lib = libs.findLibrary(libPath);

            // verify the loaded version and unload if needed
            if (lib)
            {
                // provision for manual execution of code after loading
                if (dlSymFound(lib, globalFuncName))
                {
                    loaderFunctionType function =
                        reinterpret_cast<loaderFunctionType>
                        (
                            dlSym(lib, globalFuncName)
                        );

                    if (function)
                    {
                        (*function)(true);    // force load
                    }
                    else
                    {
                        FatalIOErrorIn
                        (
                            "codedBaseClass::updateLibrary()",
                            contextDict
                        )   << "Failed looking up symbol " << globalFuncName
                            << nl << "from " << libPath << exit(FatalIOError);
                    }
                }
                else
                {
                    FatalIOErrorIn
                    (
                        "codedBaseClass::loadLibrary()",
                        contextDict
                    )   << "Failed looking up symbol " << globalFuncName << nl
                        << "from " << libPath << exit(FatalIOError);

                    lib = 0;
                    if (!libs.close(libPath, false))
                    {
                        FatalIOErrorIn
                        (
                            "codedBaseClass::loadLibrary()",
                            contextDict
                        )   << "Failed unloading library "
                            << libPath
                            << exit(FatalIOError);
                    }
                }
            }
        }
    }

    return lib;
}


void Foam::codedBaseClass::unloadLibrary
(
    const fileName& libPath,
    const string& globalFuncName,
    const dictionary& contextDict
) const
{
    void* lib = 0;

    if (libPath.empty())
    {
        return;
    }

    dlLibraryTable& libs = this->libs();

    lib = libs.findLibrary(libPath);

    if (!lib)
    {
        return;
    }

    // provision for manual execution of code before unloading
    if (dlSymFound(lib, globalFuncName))
    {
        loaderFunctionType function =
            reinterpret_cast<loaderFunctionType>
            (
                dlSym(lib, globalFuncName)
            );

        if (function)
        {
            (*function)(false);    // force unload
        }
        else
        {
            FatalIOErrorIn
            (
                "codedBaseClass::unloadLibrary()",
                contextDict
            )   << "Failed looking up symbol " << globalFuncName << nl
                << "from " << libPath << exit(FatalIOError);
        }
    }

    if (!libs.close(libPath, false))
    {
        FatalIOErrorIn
        (
            "codedBaseClass::"
            "updateLibrary()",
            contextDict
        )   << "Failed unloading library " << libPath
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::codedBaseClass::createLibrary
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    bool create = Pstream::master();

    if (create)
    {
        // Write files for new library
        if (!dynCode.upToDate(context))
        {
            // filter with this context
            dynCode.reset(context);

            this->prepare(dynCode,context);

            if (!dynCode.copyOrCreateFiles(true))
            {
                FatalIOErrorIn
                (
                    "codedBaseClass::createLibrary(..)",
                    context.dict()
                )   << "Failed writing files for" << nl
                    << dynCode.libRelPath() << nl
                    << exit(FatalIOError);
            }
        }

        if (!dynCode.wmakeLibso())
        {
            FatalIOErrorIn
            (
                "codedBaseClass::createLibrary(..)",
                context.dict()
            )   << "Failed wmake " << dynCode.libRelPath() << nl
                << exit(FatalIOError);
        }
    }


    // all processes must wait for compile to finish
    reduce(create, orOp<bool>());
}


void Foam::codedBaseClass::updateLibrary() const
{
    dynamicCode::checkSecurity
    (
        "codedBaseClass::updateLibrary()",
        dict_
    );

    dynamicCodeContext context(this->codeDict());

    // codeName: redirectType + _<sha1>
    // codeDir : redirectType
    dynamicCode dynCode
    (
        redirectType_ + context.sha1().str(true),
        redirectType_
    );
    const fileName libPath = dynCode.libPath();


    // the correct library was already loaded => we are done
    if (this->libs().findLibrary(libPath))
    {
        return;
    }

    Info<< "Using dynamicCode for " << this->description().c_str()
        << " at line " << dict_.startLineNumber()
        << " in " << dict_.name() << endl;


    // remove instantiation of fvPatchField provided by library
    this->clearRedirectPtr();

    // may need to unload old library
    unloadLibrary
    (
        oldLibPath_,
        dynamicCode::libraryBaseName(oldLibPath_),
        context.dict()
    );

    // try loading an existing library (avoid compilation when possible)
    if (!loadLibrary(libPath, dynCode.codeName(), context.dict()))
    {
        createLibrary(dynCode, context);

        loadLibrary(libPath, dynCode.codeName(), context.dict());
    }

    // retain for future reference
    oldLibPath_ = libPath;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::codedBaseClass::codedBaseClass
(
    const dictionary& dict
)
:
    dict_(dict)
{
}

Foam::codedBaseClass::codedBaseClass
()
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::codedBaseClass::~codedBaseClass()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
