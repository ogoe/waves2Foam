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

#include "stokesFifth.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(stokesFifth, 0);
addToRunTimeSelectionTable(waveTheory, stokesFifth, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

stokesFifth::stokesFifth
(
    const word& subDictName,
    const fvMesh& mesh_
)
:
    waveTheory(subDictName, mesh_),
    H_(readScalar(coeffDict_.lookup("height"))),
    h_(readScalar(coeffDict_.lookup("depth"))),
    omega_(readScalar(coeffDict_.lookup("omega"))),
    period_(2*PI_/omega_),
    phi_(readScalar(coeffDict_.lookup("phi"))),
    k_(vector(coeffDict_.lookup("waveNumber"))),
    K_(mag(k_)),

    Tsoft_(coeffDict_.lookupOrDefault<scalar>("Tsoft",period_))
{
    setCoefficients();

    checkWaveDirection(k_);
}


void stokesFifth::printCoeffs()
{
    Info << "Loading wave theory: " << typeName << endl;
}


void stokesFifth::setCoefficients()
{
    scalar S = 1.0/Foam::cosh(2*K_*h_);

    // A-coefficients
    A11_ = 1.0/sinh(K_*h_);
    A22_ = 3.0*pow(S,2.0)/(2.0*pow(1.0 - S,2.0));
    A31_ = (-4.0 - 20.0*S + 10.0*pow(S,2.0) - 13.0*pow(S,3.0))
        /(8.0*sinh(K_*h_)*pow(1.0 - S,3.0));
    A33_ = (-2.0*pow(S,2.0) + 11.0*pow(S,3.0))
        /(8.0*sinh(K_*h_)*pow(1.0 - S,3.0));
    A42_ =
        (
            12.0*S - 14.0*pow(S,2.0) - 264.0*pow(S,3.0)
          - 45.0*pow(S,4.0) - 13.0*pow(S,5.0)
        )/(24.0*pow(1 - S,5.0));
    A44_ = (10.0*pow(S,3.0) - 174.0*pow(S,4.0) + 291.0*pow(S,5.0) + 278.0*pow(S,6.0))
        /(48.0*(3.0 + 2.0*S)*pow(1 - S,5.0));
    A51_ =
        (
           - 1184.0 + 32.0*S + 13232.0*pow(S,2.0) + 21712.0*pow(S,3.0)
           + 20940.0*pow(S,4.0) + 12554.0*pow(S,5.0) - 500.0*pow(S,6.0)
           - 3341.0*pow(S,7.0) - 670.0*pow(S,8.0)
        )/(64.0*sinh(K_*h_)*(3.0 + 2.0*S)*(4.0 + S)*pow(1.0 - S,6.0));
    A53_ =
        (
            4.0*S + 105.0*pow(S,2.0) + 198.0*pow(S,3.0) - 1376.0*pow(S,4.0)
          - 1302.0*pow(S,5.0) - 117.0*pow(S,6.0) + 58.0*pow(S,7.0)
        )/(32.0*sinh(K_*h_)*(3.0 + 2.0*S)*pow(1.0 - S,6.0));
    A55_ =
        (
           - 6.0*pow(S,3.0) + 272.0*pow(S,4.0) - 1552.0*pow(S,5.0)
           + 852.0*pow(S,6.0) + 2029.0*pow(S,7.0) + 430.0*pow(S,8.0)
        )/(64.0*sinh(K_*h_)*(3.0 + 2.0*S)*(4.0 + S)*pow(1.0 - S,6.0));

    // B-coefficients
    B22_ = (1.0/tanh(K_*h_))*(1.0 + 2.0*S)/(2.0*(1.0 - S));
    B31_ = -3.0
        *(
             1.0 + 3.0*S + 3.0*pow(S,2.0)
           + 2.0*pow(S,3.0)
         )/(8.0*pow(1.0 - S,3.0));
    B42_ = (1.0/tanh(K_*h_))
        *(
             6.0 - 26.0*S - 182.0*pow(S,2.0) - 204.0*pow(S,3.0)
           - 25.0*pow(S,4.0) + 26.0*pow(S,5.0)
         )/(6.0*(3.0 + 2.0*S)*pow(1.0 - S,4.0));
    B44_ = (1.0/tanh(K_*h_))
        *(
             24.0 + 92.0*S + 122.0*pow(S,2.0) + 66.0*pow(S,3.0)
           + 67.0*pow(S,4.0) + 34.0*pow(S,5.0)
         )/(24.0*(3.0 + 2.0*S)*pow(1.0 - S,4.0));
    B53_ = 9.0
        *(
             132.0 + 17.0*S - 2216.0*pow(S,2.0) - 5897.0*pow(S,3.0)
           - 6292.0*pow(S,4.0) - 2687.0*pow(S,5.0) + 194.0*pow(S,6.0)
           + 467.0*pow(S,7.0) + 82.0*pow(S,8.0)
         )/(128.0*(3.0 + 2.0*S)*(4.0 + S)*pow(1.0 - S,6.0));
    B55_ = 5.0
        *(
             300.0 + 1579.0*S + 3176.0*pow(S,2.0) + 2949.0*pow(S,3.0)
           + 1188.0*pow(S,4.0) + 675.0*pow(S,5.0) + 1326.0*pow(S,6.0)
           + 827.0*pow(S,7.0) + 130.0*pow(S,8.0)
         )/(384.0*(3.0 + 2.0*S)*(4.0 + S)*pow(1.0 - S,6.0));

    // C-coefficients
    C0_ = sqrt(tanh(K_*h_));
    C2_ = sqrt(tanh(K_*h_))*(2.0 + 7.0*pow(S,2.0))/(4.0*pow(1.0 - S,2.0));
    C4_ = sqrt(tanh(K_*h_))
        *(
             4.0 + 32.0*S - 116.0*pow(S,2.0) - 400.0*pow(S,3.0)
           - 71.0*pow(S,4.0) + 146.0*pow(S,5.0)
         )/(32.0*pow(1.0 - S,5.0));

    // D-coefficients
    D2_ = - sqrt(1.0/tanh(K_*h_))/2.0;
    D4_ = sqrt(1.0/tanh(K_*h_))*(2.0 + 4.0*S + pow(S,2.0) + 2.0*pow(S,3.0))
        /(8.0*pow(1.0 - S,3.0));

    // E-coefficients
    E2_ = tanh(K_*h_)*(2.0 + 2.0*S + 5.0*pow(S,2.0))/(4.0*pow(1.0 - S,2.0));
    E4_ = tanh(K_*h_)
        *(
             8.0 + 12.0*S - 152.0*pow(S,2.0) - 308.0*pow(S,3.0)
           - 42.0*pow(S,4.0) + 77.0*pow(S,5.0)
         )/(32.0*pow(1.0 - S,5.0));


    // Bernoulli constant
    scalar eps(K_*H_/2.0);
    R_ = Foam::mag(g_)/K_*(0.5*Foam::sqr(C0_) + K_*h_ + Foam::sqr(eps)*E2_
    		+ Foam::pow(eps, 4.0)*E4_);

    // For K_*h_ = 0.753982 the coefficients should (approximately) fulfill:
//     Info << "A11_ = " << A11_ << " = 1.208490" << endl;
//     Info << "A22_ = " << A22_ << " = 0.799840" << endl;
//     Info << "A31_ = " << A31_ << " = -9.105340" << endl;
//     Info << "A33_ = " << A33_ << " = 0.368275" << endl;
//     Info << "A42_ = " << A42_ << " = -12.196150" << endl;
//     Info << "A44_ = " << A44_ << " = 0.058723" << endl;
//     Info << "A51_ = " << A51_ << " = 108.467921" << endl;
//     Info << "A53_ = " << A53_ << " = -6.941756" << endl;
//     Info << "A55_ = " << A55_ << " = -0.074979" << endl;
//     Info << "B22_ = " << B22_ << " = 2.502414" << endl;
//     Info << "B31_ = " << B31_ << " = -5.731666" << endl;
//     Info << "B42_ = " << B42_ << " = -32.407508" << endl;
//     Info << "B44_ = " << B44_ << " = 14.033758" << endl;
//     Info << "B53_ = " << B53_ << " = -103.445042" << endl;
//     Info << "B55_ = " << B55_ << " = 37.200027" << endl;
//     Info << "C0_  = " << C0_ << " = 0.798448" << endl;
//     Info << "C2_  = " << C2_ << " = 1.940215" << endl;
//     Info << "C4_  = " << C4_ << " = -12.970403" << endl;
//     Info << "D2_  = " << D2_ << " = -0.626215" << endl;
//     Info << "D4_  = " << D4_ << " = 3.257104" << endl;
//     Info << "E2_  = " << E2_ << " = 1.781926" << endl;
//     Info << "E4_  = " << E4_ << " = -11.573657" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalar stokesFifth::factor(const scalar& time) const
{
    scalar factor(1.0);
    if (Tsoft_ > 0.0)
    {
        factor = Foam::sin(2*PI_/(4.0*Tsoft_)*Foam::min(Tsoft_, time));
    }

    return factor;
}


scalar stokesFifth::eta
(
    const point& x,
    const scalar& time
) const
{
    scalar arg((k_ & x) - omega_*time - phi_);
    scalar eps(K_*H_/2.0);

    scalar eta = (
                    eps*Foam::cos(arg)
                  + pow(eps,2.0)*B22_*Foam::cos(2.0*arg)
                  + pow(eps,3.0)*B31_*(Foam::cos(arg) - Foam::cos(3.0*arg))
                  + pow(eps,4.0)*
                    (
                        B42_*Foam::cos(2.0*arg)
                      + B44_*Foam::cos(4.0*arg)
                    )
                  + pow(eps,5.0)*
                    (
                        - (B53_ + B55_)*Foam::cos(arg)
                        + B53_*Foam::cos(3.0*arg)
                        + B55_*Foam::cos(5.0*arg)
                    )
                 )/K_*factor(time) + seaLevel_;

    return eta;
}


scalar stokesFifth::pExcess
(
    const point& x,
    const scalar& time
) const
{
    scalar res = 0;

	// Get the reference level. Note that Eq. (29) in Fenton (1985) is
	// defined, such that z = 0 at the bottom
	scalar Z(returnZ(x) + h_);

    // Get the velocity and correct for the propagation speed (otherwise,
	// the time derivative of the velocity potential would be needed).
    vector vel = this->U(x, time);
    vel = omega_/K_*k_/K_ - vel;

    // Get the component-wise square of the velocity vector
    vector tmp = Foam::cmptMultiply(vel, vel);

    // The last component is a matter of transforming total pressure
    // into excess pressure. The equation follows from (29), Fenton (1985).
    res = R_ -Z*Foam::mag(g_) - 0.5*Foam::cmptSum(tmp)
        - ((x - referenceLevel_) & g_);
    res *= rhoWater_;

    return res;
}


//scalar stokesFifth::ddxPd
//(
//    const point& x,
//    const scalar& time,
//    const vector& unitVector
//) const
//{
//    scalar ddxPd(0);
//
//    return ddxPd;
//}


vector stokesFifth::U
(
    const point& x,
    const scalar& time
) const
{
    scalar arg((k_ & x) - omega_*time - phi_);
    scalar eps(K_*H_/2.0);
    scalar uBar( (C0_ + pow(eps,2.0)*C2_ + pow(eps,4.0)*C4_)/sqrt(K_/mag(g_)));
    scalar celerity(omega_/K_);
    scalar coeff(C0_*sqrt(mag(g_)/pow(K_,3.0)));
    scalar Z(returnZ(x) + h_);

    scalar Uhorz = celerity - uBar
        // First order
        + pow(eps,1.0)*coeff*K_*A11_*cos(arg)*cosh(K_*Z)
        // Second order
        + pow(eps,2.0)*2.0*coeff*K_*A22_*cos(2.0*arg)*cosh(2.0*K_*Z)
        // Third order
        + pow(eps,3.0)*coeff*K_*A31_*cos(arg)*cosh(K_*Z)
        + pow(eps,3.0)*3.0*coeff*K_*A33_*cos(3.0*arg)*cosh(3.0*K_*Z)
        // Fourth order
        + pow(eps,4.0)*2.0*coeff*K_*A42_*cos(2.0*arg)*cosh(2.0*K_*Z)
        + pow(eps,4.0)*4.0*coeff*K_*A44_*cos(4.0*arg)*cosh(4.0*K_*Z)
        // Fifth order
        + pow(eps,5.0)*coeff*K_*A51_*cos(arg)*cosh(K_*Z)
        + pow(eps,5.0)*3.0*coeff*K_*A53_*cos(3.0*arg)*cosh(3.0*K_*Z)
        + pow(eps,5.0)*5.0*coeff*K_*A55_*cos(5.0*arg)*cosh(5.0*K_*Z);

    Uhorz *= factor(time);

    scalar Uvert = + pow(eps,1.0)*coeff*K_*A11_*sin(arg)*sinh(K_*Z)
        // Second order
        + pow(eps,2.0)*2.0*coeff*K_*A22_*sin(2.0*arg)*sinh(2.0*K_*Z)
        // Third order
        + pow(eps,3.0)*coeff*K_*A31_*sin(arg)*sinh(K_*Z)
        + pow(eps,3.0)*3.0*coeff*K_*A33_*sin(3.0*arg)*sinh(3.0*K_*Z)
        // Fourth order
        + pow(eps,4.0)*2.0*coeff*K_*A42_*sin(2.0*arg)*sinh(2.0*K_*Z)
        + pow(eps,4.0)*4.0*coeff*K_*A44_*sin(4.0*arg)*sinh(4.0*K_*Z)
        // Fifth order
        + pow(eps,5.0)*coeff*K_*A51_*sin(arg)*sinh(K_*Z)
        + pow(eps,5.0)*3.0*coeff*K_*A53_*sin(3.0*arg)*sinh(3.0*K_*Z)
        + pow(eps,5.0)*5.0*coeff*K_*A55_*sin(5.0*arg)*sinh(5.0*K_*Z);

    Uvert *= factor(time);

    // Note "-" because of "g" working in the opposite direction
    return Uhorz*k_/K_ - Uvert*direction_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
