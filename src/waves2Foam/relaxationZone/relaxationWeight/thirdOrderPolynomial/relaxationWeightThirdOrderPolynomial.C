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

#include "relaxationWeightThirdOrderPolynomial.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationWeights
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(relaxationWeightThirdOrderPolynomial, 0);
addToRunTimeSelectionTable
(
    relaxationWeight,
    relaxationWeightThirdOrderPolynomial,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


relaxationWeightThirdOrderPolynomial::relaxationWeightThirdOrderPolynomial
(
    const word& subDictName,
    const fvMesh& mesh_
)
:
    relaxationWeight(subDictName, mesh_)
{

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void relaxationWeightThirdOrderPolynomial::computeWeights
(
    const labelList& cells,
    const scalarField& sigma,
    scalarField& weight
)
{
    forAll (weight, celli)
    {
        scalar s(1.0 - sigma[celli]);
        weight[celli] = -2.0*Foam::pow(s, 3.0) + 3.0*Foam::pow(s, 2.0);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relaxationWeights
} // End namespace Foam

// ************************************************************************* //
