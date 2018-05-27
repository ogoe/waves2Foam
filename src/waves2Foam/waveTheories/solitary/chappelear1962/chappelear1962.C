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

#include "chappelear1962.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(chappelear1962, 0);
addToRunTimeSelectionTable(waveTheory, chappelear1962, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


chappelear1962::chappelear1962
(
    const word& subDictName,
    const fvMesh& mesh_
)
:
    waveTheory(subDictName, mesh_),

    H_(readScalar(coeffDict_.lookup("height"))),
    h_(readScalar(coeffDict_.lookup("depth"))),

    propagationDirection_(vector(coeffDict_.lookup("direction"))),
    x0_(vector(coeffDict_.lookup("x0"))),

    L1_(readScalar(coeffDict_.lookup("L1"))),
    L3_(readScalar(coeffDict_.lookup("L3")))
{
    dL_ = L1_ - L3_;
    kappa_ = Foam::sqrt(3*dL_)/(2.0*h_);

    G_ = Foam::mag(g_);

    // Compute the propagation speed
    c_ = 1.0 + L3_ + dL_ + 3.0*Foam::sqr(dL_) + 5.0*dL_*L3_
        + 57.0/5.0*Foam::pow(dL_, 3.0) + 27.0*Foam::sqr(dL_)*L3_
        + 10.0*dL_*Foam::sqr(L3_);
    c_ *= Foam::sqrt(G_*h_);

    // Normalise the propagation direction
    propagationDirection_ /= Foam::mag(propagationDirection_);

    checkWaveDirection(propagationDirection_);

    // Estimate of Bernoillis constant. Appears sensitive to large amplitudes!
    point X(x0_);

    // Get surface elevation and a point on the surface
    scalar E = this->eta(X, 0.0);
    X -= (X & direction_)*direction_;
    X -= E*direction_;

    // Evaluate the velocity at the surface and correct into a steady reference
    // frame
    vector vel = this->U(X, 0.0);
    vel -= (c_*propagationDirection_);

    R_ = 0.5*Foam::cmptSum(Foam::cmptMultiply(vel, vel)) + Foam::mag(g_)*E;

}


void chappelear1962::printCoeffs()
{
    Info << "Loading wave theory: " << typeName << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar chappelear1962::shapeT
(
	const point& x,
	const scalar& time
) const
{
    scalar arg = ((x - x0_) & propagationDirection_) - c_*time;

    return Foam::tanh(kappa_*arg);
}


scalar chappelear1962::shapeS
(
	const point& x,
	const scalar& time
) const
{
    scalar arg = ((x - x0_) & propagationDirection_) - c_*time;

    return 1.0/Foam::cosh(kappa_*arg);
}


scalar chappelear1962::factor(const scalar& time) const
{
    // Dummy, as it does not make sense to ramp up a solitary wave

    return 0.0;
}


scalar chappelear1962::eta
(
    const point& x,
    const scalar& time
) const
{
	// Get the spatial shape function
	scalar T = shapeT(x, time);

	// Evaluate the function
    scalar res = 0.0;

    res += 1.0 + 2.0*dL_ + 2.0*L3_ - dL_*Foam::sqr(T) + 141.0/20.0*Foam::sqr(dL_)
        + 12.0*dL_*L3_ + Foam::sqr(L3_)
        - (5.0*Foam::sqr(dL_) + 6.0*dL_*L3_)*Foam::sqr(T)
        + 3.0/4.0*Foam::sqr(dL_)*Foam::pow(T, 4.0)
        + 8181.0/280.0*Foam::pow(dL_, 3.0) + 141.0/2.0*Foam::sqr(dL_)*L3_
        + 30.0*dL_*Foam::sqr(L3_) - Foam::sqr(T)*(2043.0/80.0*Foam::pow(dL_, 3.0)
        + 50.0*Foam::sqr(dL_)*L3_ + 15.0*dL_*Foam::sqr(L3_))
        + Foam::pow(T, 4.0)*(301.0/40.0*Foam::pow(dL_, 3.0) + 15.0/2.0*Foam::sqr(dL_)*L3_)
        - 101.0/80.0*Foam::pow(dL_, 3.0)*Foam::pow(T, 6.0);

    // Make transformation to the OF-coordinate system
    res *= h_;
    res += seaLevel_ - h_;

    return res;
}


//scalar chappelear1962::ddxPd
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


scalar chappelear1962::pExcess
(
    const point& x,
    const scalar& time
) const
{
    scalar res = 0;

/*    // Get the reference level.
    scalar Z(returnZ(x));

    // Get the velocity and correct for the propagation speed (otherwise,
    // the time derivative of the velocity potential would be needed).
    vector vel = this->U(x, time);
    vel -= (c_*propagationDirection_);

    // Get the component-wise square of the velocity vector
    vector tmp = Foam::cmptMultiply(vel, vel);

    // The last component is a matter of transforming total pressure
    // into excess pressure.

    res = R_ - Z*Foam::mag(g_) - 0.5*Foam::cmptSum(tmp)
        - (x & g_);
//    ; //- ((x - referenceLevel_) & g_);
    res *= rhoWater_;

    res -= 2.0*referencePressure();

    Info << Z << tab << R_ << tab << referencePressure() << endl;
*/
    return res;
}


vector chappelear1962::U
(
    const point& x,
    const scalar& time
) const
{
	// Get coordinate and transform to theoretical coordinate system
    scalar Z(returnZ(x));
    Z += h_;

	// Get the spatial shape functions
	scalar T = shapeT(x, time);
	scalar S = shapeS(x, time);

    // Compute the horizontal velocity component
    scalar Uhorz(0.0);

    Uhorz = 1.0 + L3_ + dL_*Foam::sqr(T) - 3.0/4.0*Foam::sqr(Z/h_)*Foam::sqr(dL_)
        *(1 - 4.0*Foam::sqr(T) + 3.0*Foam::pow(T, 4.0))
        + (2.*Foam::sqr(dL_) + 5.0*dL_*L3_)*Foam::sqr(T) + Foam::sqr(dL_)*Foam::pow(T, 4.0)
        + 3.0/16.0*Foam::pow(Z/h_, 4.0)*Foam::pow(dL_, 3.0)*
        (-2.0 + 17.0*Foam::sqr(T) - 30.0*Foam::pow(T, 4.0) + 15.0*Foam::pow(T, 6))
        - Foam::sqr(Z/h_)*(3.*Foam::pow(dL_, 3.0)/2.0 + 15.0/4.0*Foam::sqr(dL_)*L3_
        - (3.0/2.0*Foam::pow(dL_, 3.0) + 15.0*Foam::sqr(dL_)*L3_)*Foam::sqr(T)
        + (-15.0/2.0*Foam::pow(dL_, 3.0) + 45.0/4.0*Foam::sqr(dL_)*L3_)*Foam::pow(T, 4.0)
        + 15.0/2.0*Foam::pow(dL_, 3.0)*Foam::pow(T, 6.0))
        + (33.0*Foam::pow(dL_, 3.0)/5.0 + 18.0*Foam::sqr(dL_)*L3_
        + 10.0*dL_*Foam::sqr(L3_))*Foam::sqr(T)
        + (18.0/5.0*Foam::pow(dL_, 3.0) + 9.0*Foam::sqr(dL_)*L3_)*Foam::pow(T, 4.0)
        + 6.0/5.0*Foam::pow(dL_, 3.0)*Foam::pow(T, 6.0);

    // Compute the vertical velocity component
    scalar Uvert(0.0);

    Uvert = Z/(h_*h_*kappa_)*T*Foam::sqr(S)*
        (3.0*Foam::sqr(dL_)/2.0 + 3.0*Foam::sqr(Z)/(4.0*Foam::sqr(h_))*Foam::pow(dL_, 3.0)*
        (2.0 - 3.0*Foam::sqr(T)) + 3.0*Foam::pow(dL_, 3.0) + 15.0/2.0*Foam::sqr(dL_)*L3_
        + 3.0*Foam::pow(dL_, 3.0)*Foam::sqr(T));

    // Transform to dimensional units
    Uhorz *= -Foam::sqrt(G_*h_);
    Uhorz += c_;
    Uvert *= Foam::sqrt(G_*h_);

    return Uhorz*propagationDirection_ - Uvert*direction_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
