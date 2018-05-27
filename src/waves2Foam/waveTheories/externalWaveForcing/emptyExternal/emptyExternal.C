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

#include "emptyExternal.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(emptyExternal, 0);

addToRunTimeSelectionTable
(
    externalWaveForcing,
    emptyExternal,
    externalWaveForcing
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


emptyExternal::emptyExternal
(
    IOobject io,
    Time& rT,
    const fvMesh& mesh
)
:
    externalWaveForcing(io, rT, mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void emptyExternal::step()
{
    // Nothing to be done
}


scalar emptyExternal::eta
(
    const point& x,
    const scalar& time
) const
{
    return 0.0;
}


//scalar emptyExternal::ddxPd
//(
//    const point& x,
//    const scalar& time,
//    const vector& unitVector
//) const
//{
//    return 0.0;
//}


scalar emptyExternal::pExcess
(
    const point& x,
    const scalar& time
) const
{
    return 0.0;
}


vector emptyExternal::U
(
    const point& x,
    const scalar& time
) const
{
    return vector::zero;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
