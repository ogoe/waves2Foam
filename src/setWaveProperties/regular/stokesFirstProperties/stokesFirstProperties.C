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

#include "stokesFirstProperties.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(stokesFirstProperties, 0);
addToRunTimeSelectionTable(setWaveProperties, stokesFirstProperties, setWaveProperties);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

stokesFirstProperties::stokesFirstProperties
(
	const fvMesh & mesh,
	dictionary & dict,
	bool write
)
:
	setWaveProperties(mesh, dict, write)
{
	Info << "\nConstructing: " << this->type() << endl;

	period_ = readScalar( dict.lookup("period") );
	depth_  = readScalar( dict.lookup("depth") );
	omega_  = 2.0 * mathematicalConstant::pi / period_ ;
}

stokesFirstProperties::stokesFirstProperties
(
	const fvMesh & mesh,
	dictionary & dict,
	bool write,
	word string
)
:
	setWaveProperties(mesh, dict, write)
{
	Info << "\nConstructing: " << this->type() << " (Used by another wave theory)";

	period_ = readScalar( dict.lookup("period"+string) );
	depth_  = readScalar( dict.lookup("depth") );
	omega_  = 2.0 * mathematicalConstant::pi / period_;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void stokesFirstProperties::set()
{
	scalar k = linearWaveNumber();

	if ( write_ )
	{
		vector direction( vector(dict_.lookup("direction")));
		direction /= Foam::mag(direction);
		direction *= k;

		dict_.add( "waveNumber", direction, true );
		dict_.add( "omega"     , omega_   , true );
	}
}

scalar stokesFirstProperties::linearWaveNumber() const
{
	scalar lower(0.0);

    scalar upper = Foam::max( 4.0 * mathematicalConstant::pi / ( period_ * Foam::sqrt( Foam::mag(G_) * depth_)),
						      2.0 * mathematicalConstant::pi / ( Foam::pow( period_, 2.0) ) );

    scalar middle(0.5 * (lower + upper) );

    scalar valLower( Foam::pow(omega_, 2.0) - Foam::mag(G_) * lower * Foam::tanh( lower * depth_) ),
    	   valUpper( Foam::pow(omega_, 2.0) - Foam::mag(G_) * upper * Foam::tanh( upper * depth_) ),
    	   valMiddle( Foam::pow(omega_, 2.0) - Foam::mag(G_) * middle * Foam::tanh( middle * depth_) );

    while ( true )
    {
    	if ( Foam::sign( valLower ) == Foam::sign( valMiddle ) )
    	{
			lower    = middle;
			valLower = valMiddle;
    	}
    	else
    	{
    		upper    = middle;
    		valUpper = valMiddle;
    	}

    	middle = 0.5 * ( lower + upper );

    	valMiddle = Foam::pow(omega_, 2.0) - Foam::mag(G_) * middle * Foam::tanh( middle * depth_);

    	if ( Foam::mag(valMiddle) < 1.0e-13 )
    		break;
    }

	return middle;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
