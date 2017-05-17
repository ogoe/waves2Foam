/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    setWaveField

Description
    Loop over every cell in the computational domain and set VOF-ratio and
    velocity field accordingly to specified wave theory.

Author
    Niels Gjoel Jacobsen, Technical University of Denmark.  All rights reserved.

Additional information
    Implementation published and validated in the following journal article:

    @article { jacobsenFuhrmanFredsoe2011,
        Author = {Jacobsen, N G and Fuhrman, D R and Freds\o{}e, J},
        title = {{A Wave Generation Toolbox for the Open-Source CFD Library: OpenFoam\textregistered{}}},
        Journal = {{Int. J. for Numer. Meth. Fluids}},
        Year = {2012},
        Volume = {70},
        Number = {9},
        Pages = {1073-1088},
        DOI = {{10.1002/fld.2726}},
    }

\*---------------------------------------------------------------------------*/

#include "argList.H"

#if EXTBRANCH==1
    #if 310<OFVERSION
        #include "foamTime.H"
    #else
        #include "Time.H"
    #endif
#elif OFPLUSBRANCH==1
    #include "Time.H"
#else
    #include "Time.H"
#endif

//#if EXTBRANCH==1 && OFVERSION>310
//    #include "foamTime.H"
//#else
//    #include "Time.H"
//#endif

#include "fvMesh.H"
#include "volFields.H"
#include "setWaveField.H"

#include "uniformDimensionedFields.H"

#include "crossVersionCompatibility.H"
#include "externalWaveForcing.H"

using namespace Foam;

int main(int argc, char *argv[])
{

#   include "addTimeOptions.H"
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

#   include "readGravitationalAcceleration.H"

#   include "readWaveProperties.H"
#   include "createExternalWaveForcing.H"


    Info<< "\nReading field alpha\n" << endl;
    volScalarField alpha
    (
        IOobject
        (
            Foam::waves2Foam::aName(),
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field p\n" << endl;
    volScalarField pd
    (
        IOobject
        (
            Foam::waves2Foam::pName(),
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info << "Setting the wave field ...\n" << endl;

    setWaveField swf(mesh, U, alpha, pd);

    swf.correct();

    alpha.write();

    U.write();

    pd.write();

    externalWave->close();

    Info << nl << "End" << nl << endl;

    return 0;
}
