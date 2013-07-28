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

#include "lineDistribution.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(lineDistribution, 0);
addToRunTimeSelectionTable
(
    pointDistributions,
    lineDistribution,
    pointDistributions
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


lineDistribution::lineDistribution
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    pointDistributions( mesh, dict )
{
}


lineDistribution::~lineDistribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


pointField lineDistribution::evaluate()
{
    // Read needed material
    label N( readLabel(pointDict_.lookup("N")) );
    point xs( pointDict_.lookup("linestart") );
    point xe( pointDict_.lookup("lineend") );
    scalar stretch( pointDict_.lookupOrDefault<scalar>("stretch", 1.0) );

    // Define the return field
    pointField res(N, xs);

    // Compute the scaling factor
    scalar factor(0.0);

    for (int i=1; i < N; i++)
    {
        factor += Foam::pow( stretch, static_cast<scalar>(i) );
    }

    point dx( (xe - xs)/factor );

    for (int i=1; i < N; i++)
    {
        res[i] = res[i - 1] + Foam::pow( stretch, static_cast<scalar>(i) )*dx;
    }

    return res;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
