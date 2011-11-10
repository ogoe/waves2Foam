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

#include "stokesSecondProperties.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(stokesSecondProperties, 0);
addToRunTimeSelectionTable(setWaveProperties, stokesSecondProperties, setWaveProperties);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

stokesSecondProperties::stokesSecondProperties
(
	const fvMesh & mesh,
	dictionary & dict,
	bool write
)
:
	setWaveProperties(mesh, dict, write),
	sfp_( mesh, dict, false, "")
{
	Info << "\nConstructing: " << this->type() << endl;

	period_ = readScalar( dict.lookup("period") );
	depth_  = readScalar( dict.lookup("depth") );
	omega_  = 2.0 * mathematicalConstant::pi / period_ ;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void stokesSecondProperties::set()
{
	scalar k = sfp_.linearWaveNumber();

	if ( write_ )
	{
		vector direction( vector(dict_.lookup("direction")));
		direction /= Foam::mag(direction);
		direction *= k;

		dict_.add( "waveNumber", direction, true );
		dict_.add( "omega"     , omega_   , true );
	}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
