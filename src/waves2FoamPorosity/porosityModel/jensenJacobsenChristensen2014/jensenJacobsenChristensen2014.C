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

#include "jensenJacobsenChristensen2014.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(jensenJacobsenChristensen2014, 0);
addToRunTimeSelectionTable
(
    wavesPorosityModel,
    jensenJacobsenChristensen2014,
    wavesPorosityModel
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


jensenJacobsenChristensen2014::jensenJacobsenChristensen2014
(
    const fvMesh& mesh
)
:
    wavesPorosityModel(mesh),

    pZones_(mesh)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<fvMatrix<vector> > jensenJacobsenChristensen2014::ddt
(
    GeometricField<vector, fvPatchField, volMesh>& U
)
{
    return pZones_.ddt(U);
}


tmp<fvMatrix<vector> > jensenJacobsenChristensen2014::ddt
(
    const geometricOneField& rho,
    GeometricField<vector, fvPatchField, volMesh>& U
)
{
    return pZones_.ddt(rho, U);
}


tmp<fvMatrix<vector> > jensenJacobsenChristensen2014::ddt
(
    const dimensionedScalar& rho,
    GeometricField<vector, fvPatchField, volMesh>& U
)
{
	return pZones_.ddt(rho, U);
}


tmp<fvMatrix<vector> > jensenJacobsenChristensen2014::ddt
(
    const volScalarField& rho,
    GeometricField<vector, fvPatchField, volMesh>& U
)
{
    return pZones_.ddt(rho, U);
}

void jensenJacobsenChristensen2014::updatePorosity()
{
	// Store the porosity from the last time step. Needed for moving porosity
	// fields
	porosity_.storeOldTime();

	// Obtain the new porosity field from the pZones as a tmp<volScalarField>
	tmp<volScalarField> tporosity = pZones_.porosity();
	const volScalarField& poro = tporosity();

    // Set the internal field values
#if EXTBRANCH==1
	porosity_.internalField() = poro.internalField();
#elif OFPLUSBRANCH==1
    #if OFVERSION<1706
        porosity_.internalField() = poro.internalField();
    #else
        porosity_.ref() = poro.internalField();
    #endif
#else
    #if OFVERSION<400
	porosity_.internalField() = poro.internalField();
    #else
	porosity_.ref() = poro.internalField();
    #endif
#endif

	// Update boundary conditions
	porosity_.correctBoundaryConditions();
}


const volScalarField& jensenJacobsenChristensen2014::porosity() const
{
//	// Store the porosity from the last time step. Needed for moving porosity
//	// fields
////	porosity_.storeOldTime();
//
//	tmp<volScalarField> tporosity = pZones_.porosity();
//
//
//
////    return pZones_.porosity();
////    return tporosity;
    return porosity_;
}


void jensenJacobsenChristensen2014::addResistance(fvVectorMatrix& UEqn) const
{
    pZones_.addResistance(UEqn);
}


void jensenJacobsenChristensen2014::addResistance
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    pZones_.addResistance(UEqn, AU);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
