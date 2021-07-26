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

#include "waveTheory.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(waveTheory, 0);
defineRunTimeSelectionTable(waveTheory, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


waveTheory::waveTheory
(
    const word& subDictName,
    const fvMesh& mesh_
)
:
    IOdictionary
    (
        mesh_.thisDb().lookupObject<IOobject>("waveProperties")
    ),

    seaLevel_(readScalar(lookup("seaLevel"))),

#if OFPLUSBRANCH==1
    #if OFVERSION<1812
	    g_( uniformDimensionedVectorField
	        (
	            mesh_.thisDb().lookupObject<uniformDimensionedVectorField>("g")
	        ).value() ),
    #else
        g_(meshObjects::gravity::New(mesh_.thisDb().time()).value()),
    #endif
#else
    g_( uniformDimensionedVectorField
        (
            mesh_.thisDb().lookupObject<uniformDimensionedVectorField>("g")
        ).value() ),
#endif

    direction_( g_/mag(g_) ),

    coeffDict_(subDict(subDictName + "Coeffs")),

    PI_( M_PI ),

    wind_( lookupOrDefault<vector>( "wind", vector::zero ) )
{
    {
        IOdictionary transProp
        (
            IOobject
            (
                "transportProperties",
                "constant",
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        dictionary sDA(transProp.subDict(Foam::waves2Foam::airPhase()));
        rhoAir_ = (dimensionedScalar(sDA.lookup("rho"))).value();

        dictionary sDW(transProp.subDict(Foam::waves2Foam::waterPhase()));
        rhoWater_ = (dimensionedScalar(sDW.lookup("rho"))).value();
    }

    // Set the reference level
    referenceLevel_ = vector::zero;

    bool correctReference =
        this->lookupOrDefault<bool>("seaLevelAsReference", false);

    if (correctReference)
    {
        referenceLevel_ = seaLevel_*Foam::cmptMultiply(direction_, direction_);
    }
}


waveTheory::~waveTheory()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar waveTheory::referencePressure() const
{
    return referencePressure(seaLevel_);
}

scalar waveTheory::referencePressure(const scalar localSeaLevel) const
{
    return rhoWater_*Foam::mag(g_)
        *(localSeaLevel + (referenceLevel_ & direction_));
}



void waveTheory::checkWaveDirection(const vector& k) const
{
    if (Foam::mag(k & direction_) > SMALL)
    {
        FatalErrorIn("void waveTheory::checkWaveDirection(const vector& k)")
            << "The wave number " << k << " is not perpendicular to the \n"
            << "direction of the gravitational vector " << g_ << "\n"
            << endl << endl << exit(FatalError);
    }
}


void waveTheory::checkWaveDirection(const vectorField& k) const
{
    forAll (k, ki)
    {
        checkWaveDirection(k[ki]);
    }
}


scalarField waveTheory::eta
(
    const pointField& x,
    const scalar& time
) const
{
    scalarField temp(x.size(),0.0);

    forAll (x,ii)
    {
        temp[ii] = eta(x[ii],time);
    }

    return temp;
}


//scalarField waveTheory::ddxPd
//(
//    const pointField& x,
//    const scalar& time,
//    const vectorField& unitVector
//) const
//{
//    scalarField temp(x.size(),0.0);
//
//    forAll (x,ii)
//    {
//        temp[ii] = ddxPd(x[ii],time, unitVector[ii]);
//    }
//
//    return temp;
//}


vectorField waveTheory::U
(
    const pointField& x,
    const scalar& time
) const
{
    vectorField temp(x.size(),vector::zero);

    forAll (x,ii)
    {
        temp[ii] = U(x[ii],time);
    }

    return temp;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
