/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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

Class
    Foam::pressureJumpBasicFvPatchField

Description
    Foam::pressureJumpBasicFvPatchField

SourceFiles
    pressureJumpBasicFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef pressureJumpBasicFvPatchField_H
#define pressureJumpBasicFvPatchField_H

#include "jumpCyclicFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class pressureJumpBasicFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class pressureJumpBasicFvPatchField
:
    public jumpCyclicFvPatchField<Type>
{

protected:

    // Protected data

        //- "jump" field
        Field<Type> jump_;


public:


    // Constructors

        //- Construct from patch and internal field
        pressureJumpBasicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pressureJumpBasicFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given pressureJumpBasicFvPatchField onto a
        //  new patch
        pressureJumpBasicFvPatchField
        (
            const pressureJumpBasicFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        pressureJumpBasicFvPatchField
        (
            const pressureJumpBasicFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new pressureJumpBasicFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        pressureJumpBasicFvPatchField
        (
            const pressureJumpBasicFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new pressureJumpBasicFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the "jump" across the patch.
            virtual tmp<Field<Type> > jump() const
            {
                if (this->cyclicPatch().owner())
                {
                    return jump_;
                }
                else
                {
                    return refCast<const pressureJumpBasicFvPatchField<Type> >
                    (
                        this->neighbourPatchField()
                    ).jump();
                }
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "pressureJumpBasicFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
