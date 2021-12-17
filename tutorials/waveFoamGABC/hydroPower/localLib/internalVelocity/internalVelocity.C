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

#include "internalVelocity.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(internalVelocity, 0);
addToRunTimeSelectionTable(waveTheory, internalVelocity, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


internalVelocity::internalVelocity
(
    const word& subDictName,
    const fvMesh& mesh_
)
:
    waveTheory(subDictName, mesh_),

    fvMesh_(mesh_),
    boundaryName_(subDictName),

    localSeaLevel_
    (
        readScalar(coeffDict_.lookup("localSeaLevel"))
    ),

    Tsoft_(coeffDict_.lookupOrDefault<scalar>("Tsoft", 0)),

    timeIndex_(-1),

    phiBoundary_(0)
{}


void internalVelocity::printCoeffs()
{
    Info << "Loading wave theory: " << typeName << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const scalarField& internalVelocity::boundaryFlux() const
{
    if (timeIndex_ < fvMesh_.time().timeIndex())
    {
        // Update time index
        timeIndex_ = fvMesh_.time().timeIndex();

        // Get the face flux
        label patchID = fvMesh_.boundaryMesh().findPatchID(boundaryName_);
        const surfaceScalarField& phi
        (
            fvMesh_.thisDb().lookupObject<surfaceScalarField>("phi")
        );

        const volScalarField& alpha
        (
            fvMesh_.thisDb().lookupObject<volScalarField>("alpha.water")
        );

        phiBoundary_.clear();
        phiBoundary_.setSize(phi.boundaryField()[patchID].size(), 0);

        scalar Q = Foam::gSum(phi.boundaryField()[patchID]*alpha.boundaryField()[patchID]);

	const vectorField& Sf = fvMesh_.Sf().boundaryField()[patchID];
        scalar waterArea = Foam::gSum(Foam::mag(Sf)*alpha.boundaryField()[patchID]);

        phiBoundary_ = Q/waterArea*alpha.boundaryField()[patchID]*Foam::mag(Sf);
    }

    return phiBoundary_;
}


scalar internalVelocity::factor(const scalar& time) const
{
    scalar factor(1.0);

    if (Tsoft_ > 0.0)
    {
        factor = Foam::sin(2*PI_/(4.0*Tsoft_)*Foam::min(Tsoft_, time));
    }

    return factor;
}


scalar internalVelocity::eta
(
    const point& x,
    const scalar& time
) const
{
    scalar weight = factor(time);

    // Make a slow transition from the initial sea level to the boundary
    // condition
    scalar eta = weight*localSeaLevel_ + (1 - weight)*seaLevel_;

    return eta;
}


scalar internalVelocity::pExcess
(
    const point& x,
    const scalar& time
) const
{
    scalar weight = factor(time);

    scalar result = weight*rhoWater_*Foam::mag(g_)*localSeaLevel_
        + (1 - weight)*rhoWater_*Foam::mag(g_)*seaLevel_;

    result += rhoWater_*(referenceLevel_ & g_);

    return result;
}


vector internalVelocity::U
(
    const point& x,
    const scalar& time
) const
{
    // Get a reference to the velocity field
    const scalarField& phiw = boundaryFlux();

    // Get a reference to the patch ID and the face centres
    label patchID = fvMesh_.boundaryMesh().findPatchID(boundaryName_);
    const vectorField& cf = fvMesh_.Cf().boundaryField()[patchID];
    const vectorField& Sf = fvMesh_.Sf().boundaryField()[patchID];

    // Loop over all cells. Slightly inefficient, but works for now
    vector res = vector::zero;
    scalarField dist = Foam::mag(cf - x);
    scalar oldDist = GREAT;

    forAll (cf, facei)
    {
        if (dist[facei] < oldDist)
        {
            oldDist = dist[facei];

            res = phiw[facei]*Sf[facei]/Foam::magSqr(Sf[facei]);
        }
    }

    return res;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
