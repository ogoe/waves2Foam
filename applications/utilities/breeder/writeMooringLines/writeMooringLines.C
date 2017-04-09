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

Description


\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "fvMesh.H"
#include "vectorIOField.H"
#include "volFields.H"

#include "sixDoFRigidBodyMotionRestraint.H"
#include "mooringLine.H"
#include "write1DFvMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createMesh.H"

    IOdictionary dynamicMesh
    (
        IOobject
        (
            "dynamicMeshDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const dictionary& sixDict = dynamicMesh.subDict("sixDoFRigidBodyMotionCoeffs");

    // Check for restraints
    if (!sixDict.found("restraints"))
    {
        Info << "\n\nThere are no restraints specified: Ending"
             << endl << endl;

        return 0;
    }

    const dictionary& restraints = sixDict.subDict("restraints");

    // Find the mooring lines (restraints with the type 'mooringLine')
    const wordList tocs = restraints.toc();

    wordList mooringLines(tocs.size());
    label count = 0;

    forAll (tocs, toci)
    {
        if (restraints.isDict(tocs[toci]))
        {
            if
            (
                word
                (
                    restraints.subDict(tocs[toci])
                    .lookup("sixDoFRigidBodyMotionRestraint")
                ) == "mooringLine"
            )
            {
                mooringLines[count++] = tocs[toci];
            }
        }
    }

    mooringLines.setSize(count);

    if (mooringLines.size() == 0)
    {
        Info << "\n\nThere are no mooring lines specified. Ending"
             << endl << endl;
    }

    // Initialise fields for points, faces, bndFaces and cells
    labelList nPoints(mooringLines.size(), 0);
    labelList nFaces(mooringLines.size(), 0);
    labelList nBndFaces(mooringLines.size(), 0);
    labelList nCells(mooringLines.size(), 0);

    // Store reference attachment point and the anchors
    vectorField refAttachments(mooringLines.size(), vector::zero);
    vectorField anchors(mooringLines.size(), vector::zero);

    // Loop over all mooring lines
    forAll (mooringLines, resi)
    {
        // Construct the mooring line
        Foam::sixDoFRigidBodyMotionRestraints::mooringLine mL
        (
            mooringLines[resi], restraints.subDict(mooringLines[resi])
        );

        // Get reference data
        restraints.subDict(mooringLines[resi]).lookup("refAttachmentPt") >> refAttachments[resi];
        restraints.subDict(mooringLines[resi]).lookup("anchor") >> anchors[resi];

        // Set dimensions of the 1D-mesh
        mL.meshDimensions
        (
            nPoints[resi],
            nFaces[resi],
            nBndFaces[resi],
            nCells[resi]
        );
    }

       // Get the point fields from the restraints
       List<pointField> lpp(nPoints.size());

       forAll (lpp, pointi)
       {
           pointField& pp = lpp[pointi];
           pp.setSize(nPoints[pointi]);
    }

       // Initialise the displacement
    tensor initialQ = tensor::I;
    tensor Q = tensor::I;
    vector initialCentreMass(sixDict.lookup("centreOfMass"));
    vector centreMass = initialCentreMass;

    // Loop over all time folders
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << "\n" << endl;

        // Correct the centre of mass and orientation is not t = 0 s
        if (0 < runTime.timeIndex())
        {
            IOdictionary temp
            (
                IOobject
                (
                    "sixDoFRigidBodyMotionState",
                    runTime.timeName(),
                    "uniform",
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            temp.lookup("centreOfMass") >> centreMass;
            temp.lookup("orientation") >> Q;

        }

        vector totalForce = vector::zero;
        vector totalTorque = vector::zero;

        // Loop over all mooring lines
        forAll (mooringLines, resi)
        {
            // Construct the mooring line
            Foam::sixDoFRigidBodyMotionRestraints::mooringLine mL
            (
                mooringLines[resi], restraints.subDict(mooringLines[resi])
            );

            // Evaluate the instantaneous attachment
            vector attachment = centreMass
                + (Q & initialQ & (refAttachments[resi] - initialCentreMass));

            // Get the points on the 1D-mesh
            lpp[resi] = mL.centreLine(attachment, anchors[resi]);

            // Write the force at the attachment point
            vector pos(attachment); vector dummy0; vector dummy1;

            mL.restrain(pos, dummy0, dummy1);

            totalForce += dummy0;
            totalTorque += ((attachment - centreMass) ^ dummy0);
        }

        Info << "\nTotal force (mooring lines):                        " 
             << totalForce << endl
             << "Total torque in global coordinates (mooring lines): " 
             << totalTorque << endl;

        // Create the oneD fvMesh
        write1DFvMesh oneDMesh(mesh, "restraintMesh");

        // Update (create) the restraint fvMesh
        oneDMesh.updateMesh(lpp, nPoints, nFaces, nBndFaces, nCells);

        // Write the restraint fvMesh
        oneDMesh.write(true);
    }

    Info<< "\nEnd" << endl;

    return 0;
}


// ************************************************************************* //
