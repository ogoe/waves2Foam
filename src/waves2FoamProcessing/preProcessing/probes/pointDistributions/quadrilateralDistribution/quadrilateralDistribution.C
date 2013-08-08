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

#include "quadrilateralDistribution.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(quadrilateralDistribution, 0);
addToRunTimeSelectionTable
(
    pointDistributions,
    quadrilateralDistribution,
    pointDistributions
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


quadrilateralDistribution::quadrilateralDistribution
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    pointDistributions( mesh, dict )
{
}


quadrilateralDistribution::~quadrilateralDistribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


pointField quadrilateralDistribution::evaluate()
{
    // Read needed material
    label N0( readLabel(pointDict_.lookup("N0")) );
    label N1( readLabel(pointDict_.lookup("N1")) );

    point xs0( pointDict_.lookup("linestart0") );
    point xe0( pointDict_.lookup("lineend0") );
    point xe1( pointDict_.lookup("lineend1") );

    scalar stretch0( pointDict_.lookupOrDefault<scalar>("stretch0", 1.0) );
    scalar stretch1( pointDict_.lookupOrDefault<scalar>("stretch1", 1.0) );

    // Define the return field
    pointField res(N0*N1, xs0);

    // Compute the scaling factor
    scalar factor0(0.0);
    scalar factor1(0.0);

    for (int i=1; i < N0; i++)
    {
        factor0 += Foam::pow( stretch0, static_cast<scalar>(i) );
    }

    for (int i=1; i < N1; i++)
    {
        factor1 += Foam::pow( stretch1, static_cast<scalar>(i) );
    }

    point dx0( (xe0 - xs0)/factor0 );
    point dx1( (xe1 - xs0)/factor1 );

    // Compute points
    for (int j=0; j < N1; j++)
    {
        if (j != 0)
        {
            res[j*N0] = res[(j - 1)*N0]
                + Foam::pow( stretch1, static_cast<scalar>(j))*dx1;
        }

        for (int i=1; i < N0; i++)
        {
            res[i + j*N0] = res[i - 1 + j*N0]
                + Foam::pow( stretch0, static_cast<scalar>(i) )*dx0;
        }
    }

    return res;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
