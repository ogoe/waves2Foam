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

#include "cnoidalFirstProperties.H"
#include "addToRunTimeSelectionTable.H"

#include "gsl/gsl_sf_ellint.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(cnoidalFirstProperties, 0);
addToRunTimeSelectionTable
(
    setWaveProperties,
    cnoidalFirstProperties,
    setWaveProperties
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //


double lowerMBound_f ( double m, void *params )
{
    // Solves to find the value of m at which Foam::sqrt(1.0 + H/d*A) == 0.0
    struct cnoidalFirstParams *p = (struct cnoidalFirstParams * ) params;

    double d = p->depth_;
    double H = p->height_;

    double K = gsl_sf_ellint_Kcomp( Foam::sqrt(m), GSL_PREC_DOUBLE );
    double E = gsl_sf_ellint_Ecomp( Foam::sqrt(m), GSL_PREC_DOUBLE );

    double A = 2.0/m - 1.0 - 3.0/m*E/K;

    // The value of 1.0e-8 is added to ensure strictly larger than 0!
    return 1.0 - 1.0e-8 + H/d*A;
}


double cnoidalFirst_f ( double m, void *params)
{
    struct cnoidalFirstParams *p = (struct cnoidalFirstParams * ) params;

    double T = p->period_;
    double d = p->depth_;
    double G = p->g_;
    double H = p->height_;

    double K = gsl_sf_ellint_Kcomp( Foam::sqrt(m), GSL_PREC_DOUBLE );
    double E = gsl_sf_ellint_Ecomp( Foam::sqrt(m), GSL_PREC_DOUBLE );

    double A = 2.0/m - 1.0 - 3.0/m*E/K;

    return T*Foam::sqrt(G/d)*Foam::sqrt(1.0 + H/d*A)
        - Foam::sqrt( 16.0*d/(3.0*H)*m * Foam::pow(K, 2.0) );
}


double cnoidalFirstProperties::solve()
{

    int status, maxIter = 1000;
    scalar eps = 1.0e-10, m, mLower = 1.0e-15, mUpper = 1.0 - 1.0e-15;

    const gsl_root_fsolver_type *T;

    gsl_root_fsolver *s;

    gsl_function FlowerBound, F ;

    struct cnoidalFirstParams params = { d_ , H_ , T_ , G_ };

    FlowerBound.function = & lowerMBound_f;
    FlowerBound.params   = & params;

    T = gsl_root_fsolver_bisection;
    s = gsl_root_fsolver_alloc(T);

    gsl_root_fsolver_set(s, &FlowerBound, mLower, mUpper);

    for (int i = 0; i < maxIter; i++)
    {
        gsl_root_fsolver_iterate(s);
        m = gsl_root_fsolver_root(s);

        status = gsl_root_test_residual( lowerMBound_f(m, &params), eps );

        if (status == 0)
        {
            break;
        }
    }

    mLower = m;

    while (true)
    {
        if
        (
            ( cnoidalFirst_f(mLower, &params) < 0.0 &&
              cnoidalFirst_f(mUpper, &params) < 0.0 )
            ||
            ( cnoidalFirst_f(mLower, &params) > 0.0 &&
              cnoidalFirst_f(mUpper, &params)    > 0.0 )
        )
        {
            mLower = 0.999*mLower + 0.001*mUpper;
        }
        else
        {
            break;
        }

        if (Foam::mag(mLower - mUpper) < 10e-8)
        {
            return -1;
        }
    }

    F.function = &cnoidalFirst_f;
    F.params   = &params;

    gsl_root_fsolver_set(s, &F, mLower, mUpper);

    for (int i = 0; i < maxIter; i++)
    {
        gsl_root_fsolver_iterate(s);
        m = gsl_root_fsolver_root(s);

        status = gsl_root_test_residual( cnoidalFirst_f(m, &params), eps );

        if (status == 0)
        {
            break;
        }
    }

    Info << m << endl;


    return m;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


cnoidalFirstProperties::cnoidalFirstProperties
(
    const Time& rT,
    dictionary& dict,
    bool write
)
:
    setWaveProperties(rT, dict, write),
    T_( readScalar( dict.lookup("period"))),
    d_( readScalar( dict.lookup("depth"))),
    H_( readScalar( dict.lookup("height")))
{
    Info << "\nConstructing: " << this->type() << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void cnoidalFirstProperties::set( Ostream& os)
{
    scalar m = solve();

    // Write the beginning of the sub-dictionary
    writeBeginning( os );

    // Write the already given parameters
    writeGiven( os, "waveType" );

    if (dict_.found( "Tsoft" ))
    {
        writeGiven( os, "Tsoft");
    }

    writeGiven( os, "depth");
    writeGiven( os, "period");
    writeGiven( os, "height");

    if (m < 0.0)
    {
        Info << "\nPARAMETERS NOT SET\nNo cnoidal wave solution"
             << " exists for given input\n" << endl;
    }
    else
    {
        double K = gsl_sf_ellint_Kcomp( Foam::sqrt(m), GSL_PREC_DOUBLE );
        double E = gsl_sf_ellint_Ecomp( Foam::sqrt(m), GSL_PREC_DOUBLE );

        double A = 2.0/m - 1.0 - 3.0/m*E/K;

        double L =
            Foam::sqrt(16.0*m * Foam::pow(K, 2.0)*Foam::pow(d_, 3.0)/(3.0*H_));
        double c = Foam::sqrt( G_*d_*(1 + A*H_/d_));
        double omega = 2*PI_/T_;

        if (write_)
        {
            writeDerived(os, "omega", omega);
            writeDerived(os, "length", L);
            writeDerived(os, "celerity", c);

            // Locally change the write precision for m to avoid it being
            // written as 1 instead of 0.9999999999 which makes elliptic
            // integrals to infinity.
            unsigned int pre = os.precision( 14 );
            writeDerived(os, "m", m);
            os.precision( pre );
        }
    }

    writeGiven( os, "direction" );

    // Write the relaxation zone
    writeRelaxationZone( os );

    // Write the closing bracket
    writeEnding( os );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
