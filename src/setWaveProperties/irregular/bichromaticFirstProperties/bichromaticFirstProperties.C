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

#include "bichromaticFirstProperties.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(bichromaticFirstProperties, 0);
addToRunTimeSelectionTable(setWaveProperties, bichromaticFirstProperties, setWaveProperties);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

bichromaticFirstProperties::bichromaticFirstProperties
(
	const fvMesh & mesh,
	dictionary & dict,
	bool write
)
:
	setWaveProperties(mesh, dict, write),
	sfp1_( mesh, dict, write, "1"),
	sfp2_( mesh, dict, write, "2")
{
	Info << "\nConstructing: " << this->type() << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void bichromaticFirstProperties::set()
{
	scalar k1 = sfp1_.linearWaveNumber();
	scalar k2 = sfp2_.linearWaveNumber();

	if ( write_ )
	{
		vector direction1( vector(dict_.lookup("direction1")));
		vector direction2( vector(dict_.lookup("direction2")));

		direction1 /= Foam::mag(direction1);
		direction2 /= Foam::mag(direction2);

		direction1 *= k1;
		direction2 *= k2;

		dict_.add( "waveNumber1", direction1, true );
		dict_.add( "waveNumber2", direction2, true );

		dict_.add( "omega1"     , sfp1_.omega(), true );
		dict_.add( "omega2"     , sfp2_.omega(), true );
	}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
