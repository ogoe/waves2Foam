/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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

#include "overtopping.H"
#include "dictionary.H"
#include "Time.H"
#include "volFields.H"
//#include "ListListOps.H"
//#include "SortableList.H"
//#include "volPointInterpolation.H"
//#include "mapPolyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(overtopping, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        overtopping,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::overtopping::overtopping
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::regionFunctionObject(name, runTime, dict),

    runTime_(runTime),
    outputPath_(fileName::null),

    phiName_(dict.lookupOrDefault<word>("phiName","phi")),
    rhoPhiName_(dict.lookupOrDefault<word>("rhoPhiName","rho*phi")),

#if OFPLUSBRANCH==1
    #if OFVERSION > 1712
    overtoppingFilePtr_(nullptr)
    #else
    overtoppingFilePtr_(NULL)
    #endif
#else
    overtoppingFilePtr_(NULL)
#endif
{
    if (Pstream::parRun())
    {
        outputPath_ =runTime.path()/".."/"postProcessing"/name;
    }
    else
    {
        outputPath_ = runTime.path()/"postProcessing"/name;
    }

    // Remove ".."
    outputPath_.clean();

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::overtopping::~overtopping()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::overtopping::operateOnZone(const faceZone& fz) const
{
    const string& zoneName(fz.name());
    const string& ownName(this->type());
    const label N = ownName.size();

    if (!zoneName.compare(0, N, ownName))
    {
        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::overtopping::read(const dictionary& dict)
{
    // Create the output file
    if (Pstream::master())
    {
        if (overtoppingFilePtr_.empty())
        {
//            // Initialise the file
            mkDir(outputPath_ + "/" + runTime_.timeName());
            overtoppingFilePtr_.reset
            (
                new OFstream
                (
                    outputPath_ + "/" + runTime_.timeName()
                    + "/overtopping.dat"
                )
            );

            // Write header
            if (overtoppingFilePtr_.valid())
            {
                overtoppingFilePtr_() << "Time:" << tab;

                const fvMesh& mesh = runTime_.db().lookupObject<fvMesh>("region0");

                const faceZoneMesh& faceZones(mesh.faceZones());

                Info << faceZones << endl;

                forAll (faceZones, fzi)
                {
                    if (operateOnZone(faceZones[fzi]))
                    {
                        overtoppingFilePtr_() << "\t" << faceZones[fzi].name();
                    }
                }

                overtoppingFilePtr_() << endl;
            }
            else
            {
                FatalErrorIn
                (
                        "void Foam::overtopping::read( ... )"
                )
                << "Output file could not be opened in " << outputPath_
                << "/" << runTime_.timeName() << endl << endl
                << exit(FatalError);
            }
        }
    }

    return true;
}


vector Foam::overtopping::calcOvertopping
(
    const fvMesh& mesh,
    const surfaceScalarField& phi,
    const surfaceScalarField& rhoPhi,
    const faceZone& fz
) const
{
    const surfaceVectorField& Sf(mesh.Sf());
    const surfaceScalarField& magSf(mesh.magSf());

    vectorField q(fz.size(), vector::zero);

    forAll (fz, facei)
    {
        label faceI(fz[facei]);

        if (mesh.isInternalFace(faceI))
        {
            q[facei] = ((rhoPhi[faceI] - phi[faceI]*rho2_)*invRhoDiff_)
                *Sf[faceI]/magSf[faceI];
        }
        else
        {
            notImplemented("Overtopping not implemented for boundaries")
        }
    }

    vector Q(Foam::gSum(q));

    return Q;
}



bool Foam::overtopping::execute()
{
    return true;
}


bool Foam::overtopping::write()
{
    Info << "Write the overtopping" << endl;

    // Write the time with high accuracy
    if (Pstream::master())
    {
        label pres0 = overtoppingFilePtr_().precision(13);

        overtoppingFilePtr_() << runTime_.time().value() << tab;

        overtoppingFilePtr_().precision(pres0);
    }

    // Get reference to mesh and fields
    const fvMesh& mesh = runTime_.db().lookupObject<fvMesh>("region0");

    rho1_ =
        dimensionedScalar
        (
            mesh.thisDb().lookupObject<dictionary>("transportProperties")
           .subDict(Foam::waves2Foam::waterPhase()).lookup("rho")
        ).value();

    rho2_ =
        dimensionedScalar
        (
            mesh.thisDb().lookupObject<dictionary>("transportProperties")
            .subDict(Foam::waves2Foam::airPhase()).lookup("rho")
        ).value();

    invRhoDiff_ = 1.0/(rho1_ - rho2_);

    const surfaceScalarField& phi =
        mesh.lookupObject<surfaceScalarField>(phiName_);
    const surfaceScalarField& rhoPhi =
        mesh.lookupObject<surfaceScalarField>(rhoPhiName_);

    // Get reference to all face zones
    const faceZoneMesh& faceZones(mesh.faceZones());

    forAll (faceZones, fzi)
    {
        if (operateOnZone(faceZones[fzi]))
        {
            const faceZone& fZ = faceZones[fzi];

            vector Q = calcOvertopping(mesh, phi, rhoPhi, fZ);

            if (Pstream::master())
            {
                overtoppingFilePtr_() << "\t" << Q;
            }
        }
    }

    if (Pstream::master())
    {
        overtoppingFilePtr_() << endl;
    }

    return true;
}


// ************************************************************************* //
