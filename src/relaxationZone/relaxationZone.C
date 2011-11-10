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
#include "relaxationZone.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(relaxationZone, 0);

// * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * //

relaxationZone::relaxationZone
(
	const fvMesh & mesh,
	volVectorField & U,
	volScalarField & alpha
)
: 
	mesh_(mesh),
	U_(U),
	alpha_(alpha),

#if OFVERSION == 15
	relaxNames_((mesh_.db().lookupObject<IOdictionary>("waveProperties")).lookup("relaxationNames")),
#else
	relaxNames_((mesh_.thisDb().lookupObject<IOdictionary>("waveProperties")).lookup("relaxationNames")),
#endif
   	
	relaxSchemePtr_(relaxNames_.size())
{ 
	forAll (relaxNames_, relaxi)
	{
		relaxSchemePtr_[relaxi] = relaxationSchemes::relaxationScheme::New(relaxNames_[relaxi], mesh_, U_, alpha_);
	}
}

// * * * * * * * * * * * * * * * Member functions * * * * * * * * * * * * * * * * //

void relaxationZone::correct()
{	
	forAll(relaxSchemePtr_, relaxi)
	{
		relaxSchemePtr_[relaxi]->correct();
	}

	alpha_.correctBoundaryConditions();
	U_.correctBoundaryConditions();
}

tmp<volScalarField> relaxationZone::numericalBeach()
{
	// Return field
	tmp<volScalarField> tartificialViscotity
	(
		new volScalarField
		(
			IOobject
			(
				"muArti",
				mesh_.time().timeName(),
				mesh_,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh_,
			dimensionedScalar("null", dimless / dimTime, 0.0),
			"fixedValue"
		)
	);

	volScalarField & artificialViscosity( tartificialViscotity() );

	forAll(relaxSchemePtr_, relaxi)
	{
		relaxSchemePtr_[relaxi]->numericalBeach( artificialViscosity );
	}

	return tartificialViscotity;
}

} // end namespace Foam
