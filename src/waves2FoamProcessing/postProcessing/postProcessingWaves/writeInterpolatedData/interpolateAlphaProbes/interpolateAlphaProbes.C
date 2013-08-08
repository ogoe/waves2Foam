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

#include "interpolateAlphaProbes.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(interpolateAlphaProbes, 0);
addToRunTimeSelectionTable
(
    postProcessingWaves,
    interpolateAlphaProbes,
    postProcessingWaves
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


interpolateAlphaProbes::interpolateAlphaProbes
(
    const Time& rT,
    const dictionary& actionProp,
    const word& action
)
:
    postProcessingWaves( rT, actionProp, action )
{
}


interpolateAlphaProbes::~interpolateAlphaProbes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void interpolateAlphaProbes::evaluate()
{
    List<std::pair<scalar, label> > timeLabel;
    scalarField x, y, z;
    List<scalarField> alphas;

    rawAlphaProbes rap( rT_, actionProperties_, actionType_ );

    rap.readAlphaProbesData(timeLabel, x, y, z, alphas);

    scalarField t = equidistantTime(timeLabel, actionProperties_);

    scalarField weights( t.size() );
    labelList leftData( t.size() );
    labelList rightData( t.size() );

    interpolationWeights(timeLabel, t, weights, leftData, rightData);

    scalarField output( weights.size(), 0.0 );

    forAll (alphas, alphaI)
    {
        const scalarField& alpha( alphas[alphaI] );

        forAll (weights, ii)
        {
            output[ii] = weights[ii]*alpha[leftData[ii]]
                + (1.0 - weights[ii] )*alpha[rightData[ii]];
        }

        std::stringstream ss;
        ss << callName_ << "_" << alphaI;

        writeIOScalarField( output, ss.str() );
    }

    std::stringstream ss;
    ss << callName_ << "_time";

    writeIOScalarField( t, ss.str() );

    writeXYZDict( readScalar(actionProperties_.lookup("deltaT")), x, y, z);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
