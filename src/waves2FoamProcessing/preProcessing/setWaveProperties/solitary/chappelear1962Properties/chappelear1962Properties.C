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

#include "chappelear1962Properties.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(chappelear1962Properties, 0);
addToRunTimeSelectionTable
(
    setWaveProperties,
    chappelear1962Properties,
    setWaveProperties
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


chappelear1962Properties::chappelear1962Properties
(
    const Time& rT,
    dictionary& dict,
    bool write
)
:
    setWaveProperties(rT, dict, write)
{
    Info << "\nConstructing: " << this->type() << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalar chappelear1962Properties::G1
(
	const scalar L1,
	const scalar L3
) const
{
    scalar res = (L1 - L3) + 2.*L3 + 14.*Foam::sqr(L1 - L3)/5. + 6.*(L1 - L3)*L3
        + Foam::sqr(L3) + 348.*Foam::pow(L1 - L3, 3.0)/35.
        + 28.*Foam::sqr(L1 - L3)*L3 + 15.*(L1 - L3)*Foam::sqr(L3);

    return res;
}


scalar chappelear1962Properties::G2
(
	const scalar L1,
	const scalar L3,
	const scalar Hd
) const
{
    scalar res = (L1 - L3) + 17.*Foam::sqr(L1 - L3)/4. + 6.*(L1 - L3)*L3
         + 771.*Foam::pow(L1 - L3, 3.0)/40. + 85.*Foam::sqr(L1 - L3)*L3/2.
         + 15.*(L1 - L3)*Foam::sqr(L3) - Hd;

    return res;
}


scalarField chappelear1962Properties::Jacobi
(
    const scalar L1,
    const scalar L3
)
const
{
    scalarField res(4, 0.0);

    scalar dL = L1 - L3;

    res[0] = 1. + 28.*dL/5. + 6.*L3 + 3.*348./35.*Foam::sqr(dL)
        + 2.*28.*L3*dL + 15.*Foam::sqr(L3);

    res[1] = -1. + 2. - 28./5.*dL + 6.*dL - 6.*L3 + 2.*L3
        - 3.*348./35.*Foam::sqr(dL) + 28.*Foam::sqr(dL) - 56.*L3*dL
        + 2.*15.*dL*L3 - 15.*Foam::sqr(L3);

    res[2] = 1. + 2.*17./4.*dL + 6.*L3 + 3.*771./40.*Foam::sqr(dL) + 2.*85./2.*dL*L3
        + 15.*Foam::sqr(L3);

    res[3] = -1. - 17.*2./4.*dL + 6.*L1 - 12.*L3 - 3.*771./40.*Foam::sqr(dL)
        + 85./2.*Foam::sqr(dL) - 85./2.*L3*2.*dL + 15.*dL*2.*L3 - 15.*Foam::sqr(L3);

    return res;
}



void chappelear1962Properties::set(Ostream& os)
{
	// Evaluate the L1 and L3 parameters
	scalar H = readScalar(dict_.lookup("height"));
	scalar d = readScalar(dict_.lookup("depth"));

	// Provide initial conditions for the Newton-Raphson iteration
	// This follows from the lowest-order solution
	scalar L1 = H/(2*d);
	scalar L3 = -H/(2*d);
    scalar Hd = H/d;

    while (true)
    {
    	// Get the Jacobian
        scalarField J = Jacobi(L1, L3);

        // Make a right-hand side / unit matrix
        scalarField B(4.0, 0.0);
        B[0] = 1.0;
        B[3] = 1.0;

        // Perform a by-hand Gauss eliminate. Chosen to avoid cross-version
        // compatibility issues

        // First step
        scalar ratio = J[2]/J[0];
        J[2] = 0;
        J[3] -= ratio*J[1];

        B[2] -= B[0]*ratio;

        // Second step
        B[2] /= J[3];
        B[3] /= J[3];

        J[3] = 1.0;

        // Third step
        B[0] -= J[1]*B[2];
        B[1] -= J[1]*B[3];

        J[1] = 0.0;

        // Fourth step
        B[0] /= J[0];
        B[1] /= J[0];
        J[0] = 1.0;

        // Now, B contains the inverse of J

        // Update the values L1 and L3
        scalar g1 = G1(L1, L3);
        scalar g2 = G2(L1, L3, Hd);

        L1 -= B[0]*g1 + B[1]*g2;
        L3 -= B[2]*g1 + B[3]*g2;

        // Compute updated values
        g1 = G1(L1, L3);
        g2 = G2(L1, L3, Hd);

        if (Foam::mag(g1) < 1.0e-10 && Foam::mag(g2) < 1.0e-10)
        {
        	break;
        }
    }

    // Write the beginning of the sub-dictionary
    writeBeginning( os );

    // Write the already given parameters
    writeGiven(os, "waveType");

    writeGiven(os, "height");
    writeGiven(os, "depth");
    writeGiven(os, "direction");
    writeGiven(os, "x0");

    // Write the two derived parameters L1 and L3
    writeDerived(os, "L1", L1);
    writeDerived(os, "L3", L3);

    // Write the relaxation zone
    writeRelaxationZone(os);

    // Write the closing bracket
    writeEnding(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
