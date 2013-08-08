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

\*---------------------------------------------------------------------------*/

#include "interpolateForcesAndMoments.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(interpolateForcesAndMoments, 0);
addToRunTimeSelectionTable
(
    postProcessingWaves,
    interpolateForcesAndMoments,
    postProcessingWaves
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


interpolateForcesAndMoments::interpolateForcesAndMoments
(
    const Time& rT,
    const dictionary& actionProp,
    const word& action
)
:
    postProcessingWaves( rT, actionProp, action )
{
}


interpolateForcesAndMoments::~interpolateForcesAndMoments()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void interpolateForcesAndMoments::evaluate()
{
    List<std::pair<scalar, label> > timeLabel;
    vectorField forces, moments;

    rawForcesAndMoments rfm( rT_, actionProperties_, actionType_ );

    rfm.readForceAndMomentData(timeLabel, forces, moments);

    scalarField t = equidistantTime(timeLabel, actionProperties_);

    scalarField weights( t.size() );
    labelList leftData( t.size() );
    labelList rightData( t.size() );

    interpolationWeights(timeLabel, t, weights, leftData, rightData);

    vectorField output( weights.size(), vector::zero );

    // Outputs the forces as index 0
    {
        forAll (weights, ii)
        {
            output[ii] = weights[ii]*forces[leftData[ii]]
                + (1.0 - weights[ii] )*forces[rightData[ii]];
        }

        std::stringstream ss;
        ss << callName_ << "_0";

        writeIOVectorField( output, ss.str() );
    }

    // Outputs the moments as index 1
    {
        forAll (weights, ii)
        {
            output[ii] = weights[ii]*moments[leftData[ii]]
                + (1.0 - weights[ii] )*moments[rightData[ii]];
        }

        std::stringstream ss;
        ss << callName_ << "_1";

        writeIOVectorField( output, ss.str() );
    }

    // Outputs the time instances
    std::stringstream ss;
    ss << callName_ << "_time";

    writeIOScalarField( t, ss.str() );

    // Outputs the names and deltaT
    wordList wl(2);
    wl[0] = "forces";
    wl[1] = "moments";

    writeNameDict( readScalar(actionProperties_.lookup("deltaT")), wl);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
