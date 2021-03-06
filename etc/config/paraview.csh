#----------------------------------*-sh-*--------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
#    \\/     M anipulation  |
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM.
#
#     OpenFOAM is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
#     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#     for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
#
# File
#     config/paraview.csh
#
# Description
#     Setup file for paraview-3.x
#     Sourced from OpenFOAM-<VERSION>/etc/cshrc or from foamPV alias
#
# Note
#     The env. variables 'ParaView_DIR' and 'ParaView_MAJOR'
#     are required for building plugins
#------------------------------------------------------------------------------

# clean the PATH
set cleaned=`$WM_PROJECT_DIR/bin/foamCleanPath "$PATH" "$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER/cmake- $WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER/paraview-"`
if ( $status == 0 ) setenv PATH $cleaned

# determine the cmake to be used
unsetenv CMAKE_HOME
foreach cmake ( cmake-2.8.4 cmake-2.8.3 cmake-2.8.1 )
    set cmake=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER/$cmake
    if ( -r $cmake ) then
        setenv CMAKE_HOME $cmake
        setenv PATH ${CMAKE_HOME}/bin:${PATH}
        break
    endif
end

#- ParaView version, automatically determine major version:
setenv ParaView_VERSION 3.10.1
setenv ParaView_MAJOR detect


# Evaluate command-line parameters for ParaView
while ( $#argv > 0 )
    switch ($argv[1])
    case ParaView*=*:
        # name=value  -> setenv name value
        eval "setenv $argv[1]:s/=/ /"
        breaksw
    endsw
    shift
end


# set MAJOR version to correspond to VERSION
# ParaView_MAJOR is "<digits>.<digits>" from ParaView_VERSION
switch ("$ParaView_VERSION")
case "$ParaView_MAJOR".*:
    # version and major appear to correspond
    breaksw

case [0-9]*:
    # extract major from the version
    setenv ParaView_MAJOR `echo ${ParaView_VERSION} | \
        sed -e 's/^\([0-9][0-9]*\.[0-9][0-9]*\).*$/\1/'`
    breaksw
endsw


set paraviewInstDir=$WM_THIRD_PARTY_DIR/ParaView-${ParaView_VERSION}
setenv ParaView_DIR $WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER/paraview-${ParaView_VERSION}

# set paths if binaries or source are present
if ( -r $ParaView_DIR || -r $paraviewInstDir ) then
    setenv PATH ${ParaView_DIR}/bin:${PATH}
    setenv LD_LIBRARY_PATH "${ParaView_DIR}/lib/paraview-${ParaView_MAJOR}:${LD_LIBRARY_PATH}"
    setenv PV_PLUGIN_PATH $FOAM_LIBBIN/paraview-${ParaView_MAJOR}

    # add in python libraries if required
    set paraviewPython=$ParaView_DIR/Utilities/VTKPythonWrapping
    if ( -r $paraviewPython ) then
        if ($?PYTHONPATH) then
            setenv PYTHONPATH ${PYTHONPATH}:${paraviewPython}:$ParaView_DIR/lib/paraview-${ParaView_MAJOR}
        else
            setenv PYTHONPATH ${paraviewPython}:$ParaView_DIR/lib/paraview-${ParaView_MAJOR}
        endif
    endif
else
    unsetenv PV_PLUGIN_PATH
endif


unset cleaned cmake paraviewInstDir paraviewPython

# -----------------------------------------------------------------------------
