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

#include "bichromaticFirst.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(bichromaticFirst, 0);
addToRunTimeSelectionTable(waveTheory, bichromaticFirst, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


bichromaticFirst::bichromaticFirst
(
    const word& subDictName,
    const fvMesh& mesh_
)
:
    waveTheory(subDictName, mesh_),
    H1_(readScalar(coeffDict_.lookup("height1"))),
    H2_(readScalar(coeffDict_.lookup("height2"))),
    h_(readScalar(coeffDict_.lookup("depth"))),
    omega1_(readScalar(coeffDict_.lookup("omega1"))),
    omega2_(readScalar(coeffDict_.lookup("omega2"))),
    period1_(2*PI_/omega1_),
    period2_(2*PI_/omega2_),
    phi1_(readScalar(coeffDict_.lookup("phi1"))),
    phi2_(readScalar(coeffDict_.lookup("phi2"))),
    k1_(vector(coeffDict_.lookup("waveNumber1"))),
    k2_(vector(coeffDict_.lookup("waveNumber2"))),
    K1_(mag(k1_)),
    K2_(mag(k2_)),

    Tsoft_
    (
        coeffDict_.lookupOrDefault<scalar>
        (
            "Tsoft",
            Foam::max(period1_,period2_)
        )
    )
{
    checkWaveDirection(k1_);
    checkWaveDirection(k2_);
}


void bichromaticFirst::printCoeffs()
{
    Info << "Loading wave theory: " << typeName << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalar bichromaticFirst::factor(const scalar& time) const
{
    scalar factor(1.0);
    if (Tsoft_ > 0.0)
    {
        factor = Foam::sin(2*PI_/(4.0*Tsoft_)*Foam::min(Tsoft_, time));
    }

    return factor;
}


scalar bichromaticFirst::eta
(
    const point& x,
    const scalar& time
) const
{
    scalar eta = ( H1_/2.0*Foam::cos(omega1_*time - (k1_ & x) + phi1_)
                 + H2_/2.0*Foam::cos(omega2_*time - (k2_ & x) + phi2_)
                 )*factor(time)
                 + seaLevel_;
    return eta;
}


//scalar bichromaticFirst::ddxPd
//(
//    const point& x,
//    const scalar& time,
//    const vector& unitVector
//) const
//{
//
//    scalar Z(returnZ(x));
//    scalar arg1(omega1_*time - (k1_ & x) + phi1_);
//    scalar arg2(omega2_*time - (k2_ & x) + phi2_);
//
//    scalar ddxPd(0);
//
//    ddxPd =
//        (
//           rhoWater_*mag(g_)*K1_*H1_/2.0*Foam::cosh(K1_*(Z + h_))
//          /Foam::cosh(K1_*h_)*Foam::sin(arg1)
//         + rhoWater_*mag(g_)*K2_*H2_/2.0*Foam::cosh(K2_*(Z + h_))
//          /Foam::cosh(K2_*h_)*Foam::sin(arg2)
//        )*factor(time);
//
//    return ddxPd;
//}


scalar bichromaticFirst::pExcess
(
    const point& x,
    const scalar& time
) const
{
    scalar Z(returnZ(x));
    scalar arg1(omega1_*time - (k1_ & x) + phi1_);
    scalar arg2(omega2_*time - (k2_ & x) + phi2_);

    scalar res = 0;

    res += rhoWater_*mag(g_)*H1_/2.0*Foam::cosh(K1_*(Z + h_))
        /Foam::cosh(K1_*h_)*Foam::cos(arg1);
    res += rhoWater_*mag(g_)*H2_/2.0*Foam::cosh(K2_*(Z + h_))
            /Foam::cosh(K2_*h_)*Foam::cos(arg2);

    res *= factor(time);
    res += referencePressure();

    return res;
}



vector bichromaticFirst::U
(
    const point& x,
    const scalar& time
) const
{
    scalar Z(returnZ(x));

    scalar Uhorz1 = PI_*H1_/period1_ *
        Foam::cosh(K1_*(Z + h_))/Foam::sinh(K1_*h_) *
        Foam::cos(omega1_*time - (k1_ & x) + phi1_);

    scalar Uhorz2 = PI_*H2_/period2_ *
        Foam::cosh(K2_*(Z + h_))/Foam::sinh(K2_*h_) *
        Foam::cos(omega2_*time - (k2_ & x) + phi2_);

    Uhorz1 *= factor(time);
    Uhorz2 *= factor(time);

    scalar Uvert1 = - PI_*H1_/period1_ *
        Foam::sinh(K1_*(Z + h_))/Foam::sinh(K1_*h_) *
        Foam::sin(omega1_*time - (k1_ & x) + phi1_);

    scalar Uvert2 = - PI_*H2_/period2_ *
        Foam::sinh(K2_*(Z + h_))/Foam::sinh(K2_*h_) *
        Foam::sin(omega2_*time - (k2_ & x) + phi2_);

    Uvert1 *= factor(time);
    Uvert2 *= factor(time);

    // Note "-" because of "g" working in the opposite direction
    return Uhorz1*k1_/K1_ + Uhorz2*k2_/K2_ - (Uvert1 + Uvert2)*direction_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
