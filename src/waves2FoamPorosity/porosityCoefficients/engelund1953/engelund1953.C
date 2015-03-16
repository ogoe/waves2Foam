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

#include "engelund1953.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(engelund1953, 0);
addToRunTimeSelectionTable(porosityCoefficient, engelund1953, porosityCoefficient);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

engelund1953::engelund1953
(
    const dictionary & poroProp
)
:
    porosityCoefficient( poroProp )
{
    dimensionedScalar d50(poroProperties_.lookup("d50"));

    dimensionedScalar alpha(poroProperties_.lookup("alpha"));

    dimensionedScalar beta(poroProperties_.lookup("beta"));

    dimensionedScalar KC
        (
        	poroProperties_.lookupOrDefault<dimensionedScalar>
            (
        	    "KC",
        	    dimensionedScalar("KC",dimless, 10000)
        	)
        );

//    scalar poro(readScalar(poroProperties_.lookup("porosity")));
    scalar poro(readResistancePorosity(poroProperties_));

    // Compute linear resistance coefficient
    dimensionedVector d( alpha * Foam::pow3( 1 - poro ) / Foam::sqr( poro ) / Foam::sqr(d50) * vector::one);

    linearCoefficient_.value() = d.value();
    linearCoefficient_.dimensions().reset( d.dimensions() );

    // Compute quadratic resistance coefficient
    dimensionedVector f( 2.0 * beta * (1 - poro) / Foam::pow3( poro ) / d50 * vector::one);

    quadraticCoefficient_.value() = f.value();
    quadraticCoefficient_.dimensions().reset( f.dimensions() );

    Info << "Coefficients (engelund1953): " << linearCoefficient_.value() << "\t" << quadraticCoefficient_.value() << "\n" << endl;
}

engelund1953::~engelund1953()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
