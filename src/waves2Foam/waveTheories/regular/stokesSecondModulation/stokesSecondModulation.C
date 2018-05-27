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

#include "stokesSecondModulation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(stokesSecondModulation, 0);
addToRunTimeSelectionTable(waveTheory, stokesSecondModulation, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


stokesSecondModulation::stokesSecondModulation
(
    const word& subDictName,
    const fvMesh& mesh_
)
:
    waveTheory(subDictName, mesh_),
    H_(readScalar(coeffDict_.lookup("height"))),
    epsilon_(readScalar(coeffDict_.lookup("epsilon"))),
    modN_(readScalar(coeffDict_.lookup("modN"))),
    h_(readScalar(coeffDict_.lookup("depth"))),
    omega_(readScalar(coeffDict_.lookup("omega"))),
    period_(2*PI_/omega_),
    phi_(readScalar(coeffDict_.lookup("phi"))),
    k_(vector(coeffDict_.lookup("waveNumber"))),
    K_(mag(k_)),
    Tsoft_(coeffDict_.lookupOrDefault<scalar>("Tsoft",period_)),
    debug_(Switch(coeffDict_.lookup("debug")))
{
    checkWaveDirection(k_);

    if
    (
        H_/2.0 - 4.0*1.0/16.0*K_*sqr(H_)*(3.0/Foam::pow(Foam::tanh(K_*h_),3.0)
        - 1.0/Foam::tanh(K_*h_)) < 0
    )
    {
        if (debug_)
        {
            WarningIn
            (
                "label stokesSecondModulation::eta(point x, scalar time)"
            ) << endl << "The validity of stokes second order is violated."
            << endl << "a_1 < 4 a_2, being first and second order"
            << " amplitudes respectively." << endl << endl;
            Info << "a1 = " << H_/2.0 << " , a2 = "
                 << (1.0/16.0*K_*sqr(H_)*(3.0/Foam::pow(Foam::tanh(K_*h_),3.0)
                     - 1.0/Foam::tanh(K_*h_))) << endl;
        }
    }
}


void stokesSecondModulation::printCoeffs()
{
    Info << "Loading wave theory: " << typeName << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalar stokesSecondModulation::waveHeight
(
    const point& x,
    const scalar& time
) const
{
    scalar arg( ( omega_*time - (k_ & x) + phi_ )/modN_ );

    return H_*(1.0 + epsilon_*sin( arg ) );
}


scalar stokesSecondModulation::factor(const scalar& time) const
{
    scalar factor(1.0);
    if (Tsoft_ > 0.0)
    {
        factor = Foam::sin(2*PI_/(4.0*Tsoft_)*Foam::min(Tsoft_, time));
    }

    return factor;
}


scalar stokesSecondModulation::eta
(
    const point& x,
    const scalar& time
) const
{
    scalar arg(omega_*time - (k_ & x) + phi_);
    scalar modH = waveHeight( x, time );

    scalar eta = (
                     // First order contribution.
                     modH/2.0*Foam::cos(arg)
                     // Second order contribution.
                   + 1.0/16.0*K_*sqr(modH)
                   *(3.0/Foam::pow(Foam::tanh(K_*h_),3.0)
                   - 1.0/Foam::tanh(K_*h_))*Foam::cos(2.0*arg)
                 )*factor(time)  // Hot-starting.
                 + seaLevel_;      // Adding sea level.

    return eta;
}


//scalar stokesSecondModulation::ddxPd
//(
//    const point& x,
//    const scalar& time,
//    const vector& unitVector
//) const
//{
//    scalar Z(returnZ(x));
//    scalar arg(omega_*time - (k_ & x) + phi_);
//    scalar modH = waveHeight( x, time );
//
//    scalar ddxPd(0);
//
//    ddxPd = (
//                rhoWater_*mag(g_)*K_*modH/2.0*Foam::cosh(K_*(Z + h_))
//               /Foam::cosh(K_*h_)*Foam::sin(arg)
//               + 1/4*rhoWater_*mag(g_)*pow(K_,2)*pow(modH,2)
//               /Foam::sinh(2*K_*h_)
//                *( 3*Foam::cosh(2*K_*(Z + h_))/pow(Foam::sinh(K_*h_),2) - 1)
//                *Foam::sin(2*arg)
//            )*factor(time);
///*
//    Info << "ddxPd still isn't implemented. Need to think about the gradient
//    on arbitrary directed mesh with arbitrary wave number vector!
//    and arbitrary g-direction!!!" << endl;*/
//    return ddxPd;
//}


vector stokesSecondModulation::U
(
    const point& x,
    const scalar& time
) const
{
    scalar Z(returnZ(x));
    scalar cel(omega_/K_);
    scalar arg(omega_*time - (k_ & x) + phi_);
    scalar modH = waveHeight( x, time );

    // First order contribution
    scalar Uhorz = PI_*modH/period_ *
                   Foam::cosh(K_*(Z + h_))/Foam::sinh(K_*h_) *
                   Foam::cos(arg);

    // Second order contribution
    Uhorz += 3.0/16.0*cel*Foam::sqr(K_*modH)*Foam::cosh(2*K_*(Z + h_))
            /Foam::pow(Foam::sinh(K_*h_),4.0)*Foam::cos(2*arg)
             - 1.0/8.0*mag(g_)*sqr(modH)/(cel*h_);

    // First order contribution
    scalar Uvert = - PI_*modH/period_ *
                   Foam::sinh(K_*(Z + h_))/Foam::sinh(K_*h_) *
                   Foam::sin(arg);

    // Second order contribution
    Uvert += - 3.0/16.0*cel*sqr(K_*modH)*Foam::sinh(2*K_*(Z + h_))
            /Foam::pow(Foam::sinh(K_*h_), 4.0)*Foam::sin(2*arg);

    // Multiply by the time stepping factor
    Uvert *= factor(time);
    Uhorz *= factor(time);

    // Generate the velocity vector
    // Note "-" because of "g" working in the opposite direction
    return Uhorz*k_/K_ - Uvert*direction_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
