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

#include "irregular.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(irregular, 0);
addToRunTimeSelectionTable(waveTheory, irregular, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


irregular::irregular
(
    const word& subDictName,
    const fvMesh& mesh_
)
:
    waveTheory(subDictName, mesh_),
    N_(readScalar(coeffDict_.lookup("N"))),
    h_(readScalar(coeffDict_.lookup("depth"))),
    amp_("amplitude", coeffDict_, N_),
    omega_("frequency", coeffDict_, N_),
    phi_("phaselag", coeffDict_, N_),
    k_("waveNumber", coeffDict_, N_),
    K_(N_),
    compDir_(N_),
    period_(N_, 0),
    velAmp_(N_, 0),

    Tsoft_( readScalar(coeffDict_.lookup("Tsoft")))
{
    omega_ *= (2.0*PI_);

    // Compute the length of k_
    K_ = Foam::mag(k_);

    compDir_ = k_ / K_;

    // Compute the period
    forAll (period_, index)
    {
        period_[index] = 2.0*PI_/omega_[index];
    }

    // Compute the velocity amplitude
    forAll (velAmp_, index)
    {
        velAmp_[index] = PI_*2.0*amp_[index]/period_[index]
            /Foam::sinh(K_[index]*h_);
    }

    checkWaveDirection(k_);
}


void irregular::printCoeffs()
{
    Info << "Loading wave theory: " << typeName << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalar irregular::factor(const scalar& time) const
{
    scalar factor(1.0);
    if (0.0 < Tsoft_ && time < Tsoft_)
    {
        factor = Foam::sin(2*PI_/(4.0*Tsoft_)*Foam::min(Tsoft_, time));
    }

    return factor;
}


scalar irregular::eta
(
    const point& x,
    const scalar& time
) const
{
    scalar eta(0);

    forAll (amp_, index)
    {
        scalar arg = omega_[index]*time - (k_[index] & x) + phi_[index];
        eta += amp_[index]*Foam::cos(arg);
    }
    eta *= factor(time);
    eta += seaLevel_;

    return eta;
}


//scalar irregular::ddxPd
//(
//    const point& x,
//    const scalar& time,
//    const vector& unitVector
//) const
//{
//    return 0.0;
//}


scalar irregular::pExcess
(
    const point& x,
    const scalar& time
) const
{
    scalar Z(returnZ(x));

    scalar res = 0;

    forAll (amp_, index)
    {
        scalar arg0 = omega_[index]*time - (k_[index] & x) + phi_[index];
        scalar arg1 = K_[index]*(Z + h_);

        res += rhoWater_*mag(g_)*amp_[index]*Foam::cosh(arg1)
            /Foam::cosh(K_[index]*h_)*Foam::cos(arg0);
    }

    res *= factor(time);

    res += referencePressure();

    return res;
}


vector irregular::U
(
    const point& x,
    const scalar& time
) const
{
    scalar Z(returnZ(x));

    vector U(vector::zero);

    forAll (amp_, index)
    {
        scalar arg0 = omega_[index]*time - (k_[index] & x) + phi_[index];
        scalar arg1 = K_[index]*(Z + h_);

        scalar Uhorz = velAmp_[index]*Foam::cosh(arg1)*Foam::cos(arg0);
        scalar Uvert = - velAmp_[index]*Foam::sinh(arg1)*Foam::sin(arg0);

        // Note "-" because of "g" working in the opposite direction
        U += Uhorz*compDir_[index] - Uvert*direction_;
    }

    U *= factor(time);

    return U;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
