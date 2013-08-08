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

#include "circularDistribution.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(circularDistribution, 0);
addToRunTimeSelectionTable
(
    pointDistributions,
    circularDistribution,
    pointDistributions
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


circularDistribution::circularDistribution
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    pointDistributions( mesh, dict )
{
}


circularDistribution::~circularDistribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


pointField circularDistribution::evaluate()
{
    // Read needed material
    label N( readLabel(pointDict_.lookup("N")) );
    point C( pointDict_.lookup("centre") );
    scalar R( readScalar(pointDict_.lookup("radius")) );
    word axis( word( pointDict_.lookup("axis") ) );

    label c0(-1), c1(-1);

    if (axis == "x")
    {
        c0 = 1;
        c1 = 2;
    }
    else if (axis == "y")
    {
        c0 = 0;
        c1 = 2;
    }
    else if (axis == "z")
    {
        c0 = 0;
        c1 = 1;
    }
    else
    {
        FatalErrorIn("pointField circularDistribution::evaluate()")
            << "\n    The specified axis is neither x, y nor z"
            << exit(FatalError);
    }

    // Define the return field
    pointField res(N, C);

    for (int i = 0; i < N; i++)
    {
        res[i].component(c0) +=
            R*Foam::cos
            (
                (2.0*M_PI/static_cast<scalar>(N))*static_cast<scalar>(i)
            );
        res[i].component(c1) +=
            R*Foam::sin
            (
                (2.0*M_PI/static_cast<scalar>(N))*static_cast<scalar>(i)
            );
    }

    return res;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
