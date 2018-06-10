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

#include "solitaryFirst.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(solitaryFirst, 0);
addToRunTimeSelectionTable(waveTheory, solitaryFirst, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


solitaryFirst::solitaryFirst
(
    const word& subDictName,
    const fvMesh& mesh_
)
:
    waveTheory(subDictName, mesh_),
    H_(readScalar(coeffDict_.lookup("height"))),
    h_(readScalar(coeffDict_.lookup("depth"))),
    propagationDirection_(vector(coeffDict_.lookup("direction"))),
    x0_( vector(coeffDict_.lookup("x0")))
{
    G_ = Foam::mag(g_);
    c_     = Foam::sqrt( G_*(H_ + h_) );

    propagationDirection_ /= Foam::mag( propagationDirection_ );

    checkWaveDirection(propagationDirection_);
}


void solitaryFirst::printCoeffs()
{
    Info << "Loading wave theory: " << typeName << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalar solitaryFirst::factor( const scalar& time) const
{
    // Dummy, as it does not make sense to ramp up a solitary wave

    return 0.0;
}


scalar solitaryFirst::eta
(
    const point& x,
    const scalar& time
) const
{
    scalar eta = H_*1.0/Foam::pow(
        cosh(Foam::sqrt(3.0*H_/(4.0*Foam::pow(h_, 3.0)))
       *(c_*time - ((x - x0_) & propagationDirection_))) ,2.0 ) + seaLevel_;
    return eta;
}


//scalar solitaryFirst::ddxPd
//(
//    const point& x,
//    const scalar& time,
//    const vector& unitVector
//) const
//{
//    // A quite nasty expression!
//
//    return 0.0;
//}


scalar solitaryFirst::pExcess
(
    const point& x,
    const scalar& time
) const
{
    scalar Z(returnZ(x));

    scalar result( G_*H_*rhoWater_ );

    result /= (Foam::pow(Foam::cosh(Foam::sqrt(3*H_/(4.0*Foam::pow(h_, 3.0)))
        *(Foam::sqrt( G_*(h_ + H_))*time - ((x - x0_) & propagationDirection_)
        )), 4.0 ));

    result *= (2.0*(Foam::pow(h_, 3.0) + 6.0*h_*H_*Z + 3.0*H_*Foam::pow(Z, 2.0))
        + ( 2.0*Foam::pow(h_, 3.0) - 6.0*h_*H_*Z - 3.0*H_*Foam::pow(Z, 2.0))
        * Foam::cosh( Foam::sqrt( 3.0*H_/Foam::pow(h_, 3.0))*
          (Foam::sqrt(G_*(h_ + H_))*time - ((x - x0_) & propagationDirection_))
          ));

    result /= ( 4.0*Foam::pow(h_, 3.0) );

    result += referencePressure();

    return result;
}


vector solitaryFirst::U
(
    const point& x,
    const scalar& time
) const
{
    scalar Z(returnZ(x));

    scalar Uhorz(0.0);

    Uhorz = 1.0/(4.0*Foam::pow(h_, 4.0) )*Foam::sqrt(G_*h_)*H_
        *(
             2*Foam::pow(h_, 3.0) + Foam::pow(h_, 2.0)*H_
           + 12.0*h_*H_*Z + 6.0*H_*Foam::pow(Z, 2.0)
           + (
                 2*Foam::pow(h_, 3.0) - Foam::pow(h_, 2.0)*H_
               - 6.0*h_*H_*Z - 3.0*H_*Foam::pow(Z, 2.0)
             )
             *Foam::cosh
             (
                 Foam::sqrt( 3*H_/Foam::pow(h_, 3.0))
                *(
                     Foam::sqrt( G_*(h_ + H_))*time
                   - ((x - x0_) & propagationDirection_)
                 )
             )
        )
        /Foam::pow
        (
            Foam::cosh(   Foam::sqrt( 3*H_/( 4.0*Foam::pow(h_, 3.0)))
           *(
                Foam::sqrt( G_*(h_ + H_))*time
              - ((x - x0_) & propagationDirection_))
            ), 4.0
        );

    scalar Uvert(0.0);

    Uvert = 1.0/( 4.0*Foam::sqrt( G_*h_) )*Foam::sqrt(3.0)*G_
        *Foam::pow( H_/Foam::pow(h_,3.0), 1.5 )*(h_ + Z)
        *(
             2*Foam::pow(h_, 3.0) - 7.0*Foam::pow(h_, 2.0)*H_
           + 10.0*h_*H_*Z + 5.0*H_*Foam::pow(Z, 2.0)
           +(
                2*Foam::pow(h_, 3.0) + Foam::pow(h_, 2.0)*H_
              - 2.0*h_*H_*Z - H_*Foam::pow(Z, 2.0)
            )
           *Foam::cosh
            (
                Foam::sqrt( 3*H_/Foam::pow(h_, 3.0))*
                (
                    Foam::sqrt( G_*(h_ + H_))*time
                  - ((x - x0_) & propagationDirection_)
                )
            )
        )
        /Foam::pow
        (
            Foam::cosh(Foam::sqrt( 3*H_/( 4.0*Foam::pow(h_, 3.0)))
          *(
               Foam::sqrt( G_*(h_ + H_))*time
             - ((x - x0_) & propagationDirection_)
           )
        ), 4.0 )
        *Foam::tanh( Foam::sqrt( 3*H_/( 4.0*Foam::pow(h_, 3.0)))
        *(
             Foam::sqrt( G_*(h_ + H_))*time
           - ((x - x0_) & propagationDirection_)
         ));

    // Corrected sign (error) in the vertical velocity (2017.07.15, NJ)
    Uvert *= -1.0;

    return Uhorz*propagationDirection_ - Uvert*direction_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
