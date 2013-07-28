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

#include "reflectionAnalysis2DLS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(reflectionAnalysis2DLS, 0);
addToRunTimeSelectionTable
(
    postProcessingWaves,
    reflectionAnalysis2DLS,
    postProcessingWaves
);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //


void reflectionAnalysis2DLS::write
(
    const scalarField& frequencies,
    const scalarField& spectra
)
{
    Info << "        - Writing computed spectra to: " << directDir_.c_str()
         << this->type() << endl;

    mkDir( directDir_ + this->type() );

    autoPtr<OFstream> spectrumPtr_;

    // Write the right going amplitudes
    {
        std::stringstream ss;
        ss << callName_ << "_rightGoing";

        spectrumPtr_.reset
        (
            new OFstream
            (
                directDir_ + "/" + this->type() + "/"
              + ss.str() + "_spectrum.dat"
            )
        );

        for (int i = 0; i < N_; i++)
        {
            spectrumPtr_() << frequencies[2*i] << tab << spectra[4*i]
                           << endl;
            spectrumPtr_() << frequencies[2*i + 1] << tab << spectra[4*i + 1]
                           << endl;
        }
        spectrumPtr_() << frequencies[ 2*N_ ] << tab << spectra[ 4*N_ ]
                       << endl;
    }

    // Write the right going amplitudes
    {
        std::stringstream ss;
        ss << callName_ << "_leftGoing";

        spectrumPtr_.reset
        (
            new OFstream
            (
                directDir_ + "/" + this->type() + "/"
              + ss.str() + "_spectrum.dat"
            )
        );

        for (int i = 0; i < N_; i++)
        {
            spectrumPtr_() << frequencies[2*i] << tab << spectra[4*i + 2]
                           << endl;
            spectrumPtr_() << frequencies[2*i + 1] << tab << spectra[4*i + 3]
                           << endl;
        }
        spectrumPtr_() << frequencies[ 2*N_ ] << tab << spectra[ 4*N_ ]
                       << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


reflectionAnalysis2DLS::reflectionAnalysis2DLS
(
    const Time& rT,
    const dictionary& actionProp,
    const word& action
)
:
    postProcessingWaves( rT, actionProp, action ),

#   include "../../../dataDict.H"
    ,

    N_( readLabel( actionProperties_.lookup("nFreq") ) ),

    period_( readScalar( actionProperties_.lookup("period") )),

    waveNumber_( vector( actionProperties_.lookup("waveNumber")) )
{
    readIndices( dataDict_, indices_ );

    x_.setSize( indices_.size(), point::zero);

    scalarField x( dataDict_.lookup("x") );
    scalarField y( dataDict_.lookup("y") );
    scalarField z( dataDict_.lookup("z") );

    forAll (indices_ , indexi)
    {
        x_[indexi] =
            vector
            (
                x[indices_[indexi]],
                y[indices_[indexi]],
                z[indices_[indexi]]
            );
    }
}


reflectionAnalysis2DLS::~reflectionAnalysis2DLS()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void reflectionAnalysis2DLS::evaluate()
{
    if (dataType() == "scalar")
    {
        Info << "        - Computing left and right going spectra (2D)"
             << endl;

        scalar omega = 2*M_PI/period_;

        List<scalarField> input = readScalarFields( indices_ );

        scalarField time = readIOScalarField( callName_ + "_time" );

        // Creating right hand side of the over-determined system
        scalarField b( input.size()*input[0].size() );

        label count(0);

        forAll (input, inputi)
        {
            const scalarField& in( input[inputi] );

            forAll (in, ii)
            {
                b[count++] = in[ii];
            }
        }

        // Creating left hand side of the over-determined system
        List<scalarField> A( 4*N_ + 1);

        forAll (A, Ai)
        {
            scalarField& a(A[Ai]);
            a.setSize(b.size(), 1.0);
        }

        label Nin( input[0].size() );

        // Fill the A-matrix
        for (int nf=0; nf < N_; nf++)
        {
            scalar J( nf + 1.0 );
            count = 0;

            scalarField& a0( A[ 4*nf ] );
            scalarField& a1( A[ 4*nf + 1 ] );
            scalarField& a2( A[ 4*nf + 2 ] );
            scalarField& a3( A[ 4*nf + 3 ] );

            forAll (x_, xi)
            {
                point X( x_[xi] );

                for (int n=0; n<Nin; n++)
                {
                    a0[count] =
                        Foam::cos(J*(omega*time[n] - (waveNumber_ & X)));
                    a1[count] =
                        Foam::sin(J*(omega*time[n] - (waveNumber_ & X)));
                    a2[count] =
                        Foam::cos(J*(omega*time[n] + (waveNumber_ & X)));
                    a3[count] =
                        Foam::sin(J*(omega*time[n] + (waveNumber_ & X)));

                    count++;
                }
            }
        }

        spectralMethodsLeastSquaresBased smls( rT_, actionProperties_ );

        smls.solve( A, b );

        scalarField frequencies( smls.frequencies( N_ ) );

        write( frequencies, b);
    }
    else
    {
        FatalErrorIn
        (
            "void reflectionAnalysis2DFFT::evaluate()"
        )   << "Only supports input of scalar quantities (surface elevation) "
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
