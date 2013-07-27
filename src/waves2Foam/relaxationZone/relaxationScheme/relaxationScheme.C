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

#include "relaxationScheme.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationSchemes
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(relaxationScheme, 0);
defineRunTimeSelectionTable(relaxationScheme, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

relaxationScheme::relaxationScheme
(
    const word& subDictName,
    const fvMesh& mesh,
    vectorField& U,
    scalarField& alpha
)
:
    IOdictionary
    (
        mesh.thisDb().lookupObject<IOobject>("waveProperties")
    ),
    convexPolyhedral(mesh, true),
    mesh_(mesh),
    U_(U),
    alpha_(alpha),
    coeffDict_(subDict(subDictName + "Coeffs").subDict("relaxationZone"))
{
    relaxShape_  = relaxationShapes::relaxationShape::New(subDictName, mesh_);
    relaxWeight_ = relaxationWeights::relaxationWeight::
        New(subDictName, mesh_);
    waveProps_   = waveTheories::waveTheory::New(subDictName, mesh_);
    numBeach_    = numericalBeaches::numericalBeach::New(subDictName, mesh_ );
}


relaxationScheme::~relaxationScheme()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void relaxationScheme::signedPointToSurfaceDistance
(
    const pointField& pp,
    scalarField& sd
)
{
    forAll (pp, pointi)
    {
        sd[pointi] = signedPointToSurfaceDistance(pp[pointi]);
    }
}


scalar relaxationScheme::signedPointToSurfaceDistance
(
    const point& pp
) const
{
    scalar temp = waveProps_->eta(pp, db().time().value() );
    temp += ( waveProps_->returnDir() & pp );
    temp *= -1.0;

    return temp;
}


void relaxationScheme::numericalBeach
(
    volScalarField& artVisc
)
{
    const labelList& cc( cells() );
    const scalarField& ss( sigma() );

    numBeach_->correct( cc, ss, artVisc );

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relaxationSchemes
} // End namespace Foam

// ************************************************************************* //
