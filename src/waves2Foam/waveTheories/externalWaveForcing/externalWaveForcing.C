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

#include "externalWaveForcing.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(externalWaveForcing, 0);
defineRunTimeSelectionTable(externalWaveForcing, externalWaveForcing);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


externalWaveForcing::externalWaveForcing
(
    IOobject io,
    Time& rT,
    const fvMesh& mesh
)
:
    regIOobject(io),

    rT_(rT),

    mesh_(mesh)
{

}


externalWaveForcing::~externalWaveForcing()
{}


autoPtr<externalWaveForcing> externalWaveForcing::New
(
    IOobject io,
    Time& rT,
    const fvMesh& mesh
)
{
    word externalType
    (
        (io.db().lookupObject<IOdictionary>("waveProperties"))
        .lookupOrDefault<word>("externalForcing", "emptyExternal")
    );

#if OFVERSION < 2206
    externalWaveForcingConstructorTable::iterator cstrIter =
        externalWaveForcingConstructorTablePtr_->find(externalType);

    if (cstrIter == externalWaveForcingConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "externalWaveForcing::New(IOobject, const Time&)"
        )   << "Unknown type of external wave forcing: " << externalType
            << endl << endl
            << "Valid external forcings are :" << endl
            << externalWaveForcingConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<externalWaveForcing>(cstrIter()(io, rT, mesh));
#else
   auto* cstrIter = externalWaveForcingConstructorTable(externalType);

    if (!cstrIter)
    {
        FatalErrorIn
        (
                "externalWaveForcing::New(IOobject, const Time&)"
        )   << "Unknown type of external wave forcing: " << externalType
        << endl << endl
        << "Valid external forcings are :" << endl
        << externalWaveForcingConstructorTablePtr_->toc()
        << exit(FatalError);
    }

    return autoPtr<externalWaveForcing>(cstrIter(io, rT, mesh));
#endif
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //





// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
