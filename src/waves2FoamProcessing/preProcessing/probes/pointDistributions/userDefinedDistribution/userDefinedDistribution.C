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

#include "userDefinedDistribution.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(userDefinedDistribution, 0);
addToRunTimeSelectionTable
(
    pointDistributions,
    userDefinedDistribution,
    pointDistributions
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


userDefinedDistribution::userDefinedDistribution
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    pointDistributions( mesh, dict )
{
}


userDefinedDistribution::~userDefinedDistribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


pointField userDefinedDistribution::evaluate()
{
    // Read needed material
    scalarField x("xValues", pointDict_, readLabel( pointDict_.lookup("N")));
    scalarField y("yValues", pointDict_, readLabel( pointDict_.lookup("N")));
    scalarField z("zValues", pointDict_, readLabel( pointDict_.lookup("N")));

    // Define the return field
    pointField res(x.size(), point::zero);

    forAll (res, pointi)
    {
        res[pointi] = point(x[pointi], y[pointi], z[pointi]);
    }

    return res;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
