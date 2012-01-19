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

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

waveTheory::waveTheory
(
	const word & subDictName,
	const fvMesh & mesh_
)
:
    IOdictionary
    (
#if OFVERSION == 15
        mesh_.db().lookupObject<IOobject>("waveProperties")
#else
        mesh_.thisDb().lookupObject<IOobject>("waveProperties")
#endif
    ),

    seaLevel_(readScalar(lookup("seaLevel"))),

// Takes care of the fact that the gravity vector is defined differently between OF1.5 and OF1.6+
#if OFVERSION == 15
    g_( dimensionedVector( mesh_.db().lookupObject<IOdictionary>("environmentalProperties").lookup("g") ).value() ),
#else
	g_( uniformDimensionedVectorField( mesh_.thisDb().lookupObject<uniformDimensionedVectorField>("g")).value() ),
#endif
    direction_( g_ / mag(g_) ),

    coeffDict_(subDict(subDictName + "Coeffs")),

    PI_( 4.0 * Foam::atan(1.0) )
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
        dictionary sD(transProp.subDict("phase1"));
        rhoWater_ = (dimensionedScalar(sD.lookup("rho"))).value();
    }
}


waveTheory::~waveTheory()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalarField waveTheory::eta
(
	const pointField & x,
	const scalar & time
)
{
    scalarField temp(x.size(),0.0);

    forAll(x,ii)
    {
        temp[ii] = eta(x[ii],time);
    }

    return temp;
}

scalarField waveTheory::ddxPd
(
	const pointField & x,
	const scalar & time,
	const vectorField & unitVector
)
{
    scalarField temp(x.size(),0.0);

    forAll(x,ii)
    {
        temp[ii] = ddxPd(x[ii],time, unitVector[ii]);
    }

    return temp;
}

vectorField waveTheory::U
(
	const pointField & x,
	const scalar & time
)
{
    vectorField temp(x.size(),vector::zero);

    forAll(x,ii)
    {
        temp[ii] = U(x[ii],time);
    }

    return temp;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
