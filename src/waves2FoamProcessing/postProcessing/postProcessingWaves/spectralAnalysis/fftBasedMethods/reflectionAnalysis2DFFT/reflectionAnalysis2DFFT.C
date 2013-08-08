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

#include "reflectionAnalysis2DFFT.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(reflectionAnalysis2DFFT, 0);
addToRunTimeSelectionTable
(
    postProcessingWaves,
    reflectionAnalysis2DFFT,
    postProcessingWaves
);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //


scalarField reflectionAnalysis2DFFT::linearWaveNumbers
(
    const scalarField& frequencies
)
{
    // Prepare return field
    scalarField res( frequencies.size(), 0.0);

    // Load the stokes first calculator
    stokesFirstProperties stp( rT_ , dataDict_ );

    // Compute the linear wave number for each freqency
    forAll (frequencies, freqi)
    {
        res[freqi] = stp.linearWaveNumber( depth_, frequencies[freqi] );
    }

    return res;
}


void reflectionAnalysis2DFFT::writeReflectionIncident
(
    const scalarField& frequencies,
    const scalarField& spectrumRight,
    const scalarField& spectrumLeft,
    const scalarField& determinant
)
{
    Info << "        - Writing computed spectra to: " << directDir_.c_str()
         << this->type() << endl;

    mkDir( directDir_ + this->type() );

    autoPtr<OFstream> spectrumPtr_;

    // Write left going spectra
    {
        std::stringstream ss;
        ss << callName_ << "_" << writeIndex_ << "_leftGoing";

        spectrumPtr_.reset
        (
            new OFstream
            (
                directDir_ + "/" + this->type() + "/"
              + ss.str() + "_spectrum.dat"
            )
        );

        forAll (frequencies, freqi)
        {
            spectrumPtr_() << frequencies[freqi] << tab << spectrumLeft[freqi]
                           << endl;
        }
    }

    // Write right going spectra
    {
        std::stringstream ss;
        ss << callName_ << "_" << writeIndex_ << "_rightGoing";

        spectrumPtr_.reset
        (
            new OFstream
            (
                directDir_ + "/" + this->type() + "/"
              + ss.str() + "_spectrum.dat"
            )
        );

        forAll (frequencies, freqi)
        {
            spectrumPtr_() << frequencies[freqi] << tab << spectrumRight[freqi]
                           << endl;
        }
    }

    // Write determinant
    {
        std::stringstream ss;
        ss << callName_ << "_" << writeIndex_ << "_determinant";

        spectrumPtr_.reset
        (
            new OFstream
            (
                directDir_ + "/" + this->type() + "/" + ss.str()
              + "_spectrum.dat"
            )
        );

        forAll (frequencies, freqi)
        {
            spectrumPtr_() << frequencies[freqi] << tab << determinant[freqi]
                           << endl;
        }
    }
}


void reflectionAnalysis2DFFT::decomposeAmplitudes
(
    const scalarField& k,
    const List<Field<complex> >& amps,
    Field<complex>& ampRight,
    Field<complex>& ampLeft,
    scalarField& determinant
)
{
    // Number of wave gauges
    label P( X_.size() );

    // Complex 'i'.
    complex ii( 0, 1 );

    // Set the correct size of return fields
    if (ampRight.size() != k.size())
    {
        ampRight.setSize(k.size(), complex::zero);
    }

    if (ampLeft.size() != k.size())
    {
        ampLeft.setSize(k.size(), complex::zero);
    }

    if (determinant.size() != k.size())
    {
        determinant.setSize(k.size(), 0.0);
    }

    // Loop over all frequencies (wave numbers)
    forAll (k, freqi)
    {
        scalar kj( k[freqi] );

        // Compute the weight (heuristic suggestion in article)
        // (Eq. 21 and eq. 22)
        scalarField Wj( X_.size(), 0.0);

        for (int p = 0; p < P; p++)
        {
            for (int q = 0; q < P; q++)
            {
                scalar delta( kj*( X_[p] - X_[q] ));
                Wj[p] += Foam::sqr(Foam::sin(delta))
                    /(1 + Foam::sqr(delta/M_PI));
            }
        }

        // Compute the denominator D (Eq. 12c)
        scalar D(0);
        for (int p=0; p<P; p++)
        {
            for (int q=0; q < p; q++)
            {
                D +=
                    (
                        4.0*Wj[p]*Wj[q]
                       *Foam::sqr( Foam::sin(kj*(X_[p] - X_[q])))
                    );
            }
        }

        determinant[freqi] = D;

        // Compute the complex coefficients C (Eq. 15)
        Field<complex> C( X_.size(), complex::zero );

        for (int p=0; p < P; p++)
        {
            C[p] = complex::zero;

            for (int q=0; q < P; q++)
            {
                C[p] +=
                     (
                         Wj[q]*Foam::sin( kj*( X_[p] - X_[q] ))
                        *exp( ii*kj*(X_[q] - X_[0]))
                     );
            }

            C[p] *= (2.0*ii*Wj[p]*exp( ii*kj*X_[0] )/D);
        }

        // Compute the left and right going discrete fourier coefficients.
        // (Eq. 14a and 14b)
        for (int p = 0; p < P; p++)
        {
            const Field<complex>& amp( amps[p] );

            ampLeft[freqi]  += ( (C[p].conjugate())*amp[freqi] );
            ampRight[freqi] += ( C[p]*amp[freqi] );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


reflectionAnalysis2DFFT::reflectionAnalysis2DFFT
(
    const Time& rT,
    const dictionary& actionProp,
    const word& action
)
:
    postProcessingWaves( rT, actionProp, action ),

#   include "../../../dataDict.H"
{
    readIndices( dataDict_, indices_ );

    deltaT_ = readDeltaT( dataDict_ );

    depth_ = readScalar( actionProperties_.lookup("depth") );

    coordName_ = word(actionProperties_.lookup("coordName"));

    writeIndex_ = actionProperties_.lookupOrDefault<label>("writeIndex", 0);

    X_.setSize( indices_.size() );

    scalarField allX = dataDict_.lookup( coordName_ );

    forAll (indices_, indexi)
    {
        X_[indexi] = allX[indices_[indexi]];
    }
}


reflectionAnalysis2DFFT::~reflectionAnalysis2DFFT()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void reflectionAnalysis2DFFT::evaluate()
{
    if (dataType() == "scalar")
    {
        Info << "        - Computing left and right going spectra (2D)"
             << endl;

        spectralMethodsFFTBased smfft(rT_, actionProperties_);

        List<scalarField> input = readScalarFields(indices_);

        smfft.initSweep(input[0]);

        scalarField frequencies = smfft.frequencies(deltaT_);

        scalarField k = linearWaveNumbers(frequencies);

        scalarField specLeft(k.size(), 0.0), specRight(k.size(), 0.0);
        scalarField determinant(k.size(), 0.0);

        scalar factor = 2.0*deltaT_
            /static_cast<scalar>( smfft.nSweeps()*smfft.nBins() );

        while (true)
        {
            // Make the Fourier transform of the next bit of the data set
            List<Field<complex> > transforms = smfft.fft( input );

            // In case the Fourier transform sequence is done, a List of size
            // zero is returned.
            if (!transforms.size())
            {
                break;
            }

            // Prepare decomposition FFT fields
            Field<complex> ampRight;
            Field<complex> ampLeft;

            // Decompose amplitudes
            decomposeAmplitudes(k, transforms, ampRight, ampLeft, determinant);

            // Add the decomposed DFTs to the left and right going power
            // spectra
            forAll (ampRight, freqi)
            {
                specRight[freqi] +=
                    (
                        factor
                       *(ampRight[freqi]*ampRight[freqi].conjugate()).Re()
                    );
                specLeft[freqi]  +=
                    (factor*(ampLeft[freqi]*ampLeft[freqi].conjugate()).Re());
            }
        }

        writeReflectionIncident(frequencies, specRight, specLeft, determinant);
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
