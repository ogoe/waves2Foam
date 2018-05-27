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

#include "combinedWaves.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(combinedWaves, 0);
addToRunTimeSelectionTable(waveTheory, combinedWaves, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


combinedWaves::combinedWaves
(
    const word& subDictName,
    const fvMesh& mesh_
)
:
    waveTheory(subDictName, mesh_),

    combinedWavesNames_(coeffDict_.lookup("combinedWaveNames")),
    combinedWavesPtr_(combinedWavesNames_.size())
{
    if (combinedWavesNames_.size() == 0)
    {
        FatalErrorIn
            (
             "Foam::waveTheories::combinedWaves(const word& subDictName, const fvMesh& mesh_)"
            )   << "The size of the combining waves is "
            << combinedWavesNames_.size() << endl << endl
            << "There should be at least one (1) wave type." << endl
            << exit(FatalError);
    }

    forAll (combinedWavesPtr_, cI)
    {
        combinedWavesPtr_[cI] = waveTheories::waveTheory::
            New(combinedWavesNames_[cI], mesh_);
    }
}


void combinedWaves::printCoeffs()
{
    Info << "Loading wave theory: " << typeName << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar combinedWaves::factor(const scalar& time) const
{
    // Not used in the present case, as the factor is multiplied on to the data
    // in the individual wave theories.
    // Needed to be here, as it is an abstract part of waveTheories::waveTheory
    return 1;
}


scalar combinedWaves::eta
(
    const point& x,
    const scalar& time
) const
{
    scalar eta(0);

    forAll (combinedWavesPtr_, cI)
    {
        eta += combinedWavesPtr_[cI]->eta(x, time);
    }

    // There must only be corrected for the seaLevel_ once.
    eta -= ( (combinedWavesNames_.size() - 1.0)*seaLevel_ );

    return eta;
}


//scalar combinedWaves::ddxPd
//(
//    const point& x,
//    const scalar& time,
//    const vector& unitVector
//) const
//{
//    scalar ddxPd(0.0);
//
//    forAll (combinedWavesPtr_, cI)
//    {
//        ddxPd += combinedWavesPtr_[cI]->ddxPd(x, time, unitVector);
//    }
//
//    return ddxPd;
//}


scalar combinedWaves::pExcess
(
    const point& x,
    const scalar& time
) const
{
    scalar res(0);

    forAll (combinedWavesPtr_, cI)
    {
        res += combinedWavesPtr_[cI]->pExcess(x, time);
    }

    // There must only be corrected for the reference pressure once.
    res -= ((combinedWavesNames_.size() - 1.0)*referencePressure());

    return res;
}


vector combinedWaves::U
(
    const point& x,
    const scalar& time
) const
{
    vector U(vector::zero);

    forAll (combinedWavesPtr_, cI)
    {
        U += combinedWavesPtr_[cI]->U(x, time);
    }


    return U;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
