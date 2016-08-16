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
    relaxationZoneLayout

Description
    Loop over every cell in the computational domain and set the relaxation zone
    coordinate sigma and show which cells belong to which relaxatio zone.

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.

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


#include "fvMesh.H"
#include "volFields.H"
#include "relaxationShape.H"
#include "relaxationWeight.H"

#include "uniformDimensionedFields.H"

using namespace Foam;

int main(int argc, char *argv[])
{

#   include "addTimeOptions.H"
#   include "setRootCase.H"

#   include "createTime.H"

#   include "createMesh.H"

#   include "readGravitationalAcceleration.H"

#   include "readWaveProperties.H"

    Info<< "Creating field relaxationZoneLayout\n" << endl;

    volScalarField rzl
    (
        IOobject
        (
            "relaxationZoneLayout",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("null", dimless, -1),
        "zeroGradient"
    );

    Info<< "Creating field relaxationZoneLayoutSigma\n" << endl;

    volScalarField srzl
    (
        IOobject
        (
            "relaxationZoneSigmaValue",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("null", dimless, 0),
        "zeroGradient"
    );

    volScalarField wrzl
    (
        IOobject
        (
            "relaxationZoneWeightOnComputed",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("null", dimless, 1.0),
        "zeroGradient"
    );

    // Read the relaxation names from waveProperties
    wordList relaxNames(waveProperties.lookup("relaxationNames"));

    // Create the relaxation shape function and add the data to the above
    // generated fields
    forAll (relaxNames, relaxi)
    {
        // Get relaxation zone cells and coorsponding sigma coordinates
        autoPtr<relaxationShapes::relaxationShape> relaxShape
            = relaxationShapes::relaxationShape::New(relaxNames[relaxi], mesh);

        const labelList& cells(relaxShape->cells());

        const scalarField& sigma(relaxShape->sigma());

        forAll (cells, celli)
        {
            rzl[cells[celli]] += (relaxi + 1);
            srzl[cells[celli]] = sigma[celli];
        }

        Info << relaxNames[relaxi] << " has "
             << cells.size() << " cells\n" << endl;

        // Compute relaxation zone weights
        autoPtr<relaxationWeights::relaxationWeight> relaxWeight
            = relaxationWeights::relaxationWeight::
              New(relaxNames[relaxi], mesh);

        scalarField weights(sigma.size(), 1.0);

        relaxWeight->weights(cells, sigma, weights);

        forAll (weights, celli)
        {
            wrzl[cells[celli]] = weights[celli];
        }
    }

    Info << "Write the fields to Time = " << runTime.timeName() << endl;
    rzl.write();
    srzl.write();
    wrzl.write();

    Info << nl << "End" << endl;

    return 0;
}
