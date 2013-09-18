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

#include "spectralMethodsFFTBased.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(spectralMethodsFFTBased, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //


void Foam::spectralMethodsFFTBased::checkBins()
{
    if ((bins_ % 2) == 1)
    {
        FatalErrorIn("void Foam::spectralMethodsFFTBased::checkBins()")
                << "The number of frequency bins (" << bins_
                << ") given in the dictionary" << endl
                << "is not an even number" << endl << exit(FatalError);
    }
}


Field<complex> spectralMethodsFFTBased::computeFourierTransform
(
    const scalarField& input
)
{
    Field<complex> res( bins_/2 );

    // Prepare data for sweep i
    double data[bins_];

    for (int m = 0; m < bins_; m++)
    {
        data[m] = input[ sweepCount_*step_ + m ];
    }

    // Compute the FFT
    gsl_fft_real_transform (data, 1, bins_, real, work);

    for (int m = 1; m < bins_/2; m++ )
    {
        res[m - 1].Re() = data[2*m - 1];
        res[m - 1].Im() = data[2*m];
    }

    res[bins_/2 - 1] = data[ bins_ - 1];

    return res;
}


void Foam::spectralMethodsFFTBased::powerSpectrum
(
    const scalarField& input,
    const scalar& deltaT,
    scalarField& spectrum
)
{
    // Getting the sweep parameters
    resetSweep();
    initSweep( input );

    // GSL-FFT data
    scalar factor( 2.0*deltaT/static_cast<scalar>( sweeps_*bins_ ) );

    for (sweepCount_ = 0; sweepCount_ < sweeps_; sweepCount_++)
    {
        Field<complex> transform( computeFourierTransform( input ) );

        for (int m=0; m < bins_/2; m++)
        {
            spectrum[m] += (factor*transform[m]*transform[m].conjugate()).Re();
        }
    }
}


void spectralMethodsFFTBased::powerSpectrum
(
    const vectorField& input,
    const scalar& deltaT,
    vectorField& spectrum
)
{
    // Compute the spectrum for each vector component
    for (int i=0; i<3; i++)
    {
        scalarField comp = input.component(i);
        scalarField speci = spectrum.component(i);

        powerSpectrum(comp, deltaT, speci);

        spectrum.replace( i, speci);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


spectralMethodsFFTBased::spectralMethodsFFTBased
(
    const Time& mesh,
    const dictionary& actionProp
)
:
    dict_( actionProp ),

    bins_( readLabel( dict_.lookup("freqBins") ) ),

    sweeps_( -1 ),

    sweepCount_( -1 )
{
    checkBins();

    scalar overlap( readScalar( dict_.lookup("windowShiftFraction") ) );

    step_ = static_cast<label>( bins_*overlap );

    if (overlap <= 0 || step_ == 0)
    {
        FatalErrorIn
        (
            "void Foam::spectralMethodsFFTBased::spectralMethodsFFTBased(const fvMesh&, const dictionary&)"
        )
        << "The overlap-factor is negative or the overlap is so small"
        << " that the resulting step" << endl
        << "in window size is 0" << endl << exit(FatalError);
    }

    // Allocation memory for FFT
    real = gsl_fft_real_wavetable_alloc(bins_);
    work = gsl_fft_real_workspace_alloc(bins_);
}


spectralMethodsFFTBased::~spectralMethodsFFTBased()
{
    // Free GSL-FFT
    gsl_fft_real_wavetable_free (real);
    gsl_fft_real_workspace_free (work);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


List<scalarField> spectralMethodsFFTBased::powerSpectra
(
    const List<scalarField>& inputData,
    const scalar& deltaT
)
{
    List<scalarField> spectra( inputData.size() );

    forAll (spectra, speci)
    {
        scalarField& spectrum( spectra[speci] );
        spectrum.setSize( bins_/2, 0.0 );

        const scalarField& input( inputData[ speci ] );

        powerSpectrum( input, deltaT, spectrum );
    }

    Info << "          (The power spectra is averaged over " << sweeps_ << " individual spectra)." << endl;

    return spectra;
}


List<vectorField> spectralMethodsFFTBased::powerSpectra
(
    const List<vectorField>& inputData,
    const scalar& deltaT
)
{
    List<vectorField> spectra( inputData.size() );

    forAll (spectra, speci)
    {
        vectorField& spectrum( spectra[speci] );
        spectrum.setSize( bins_/2, vector::zero );

        const vectorField& input( inputData[ speci ] );

        powerSpectrum( input, deltaT, spectrum );
    }

    Info << "          (The power spectra is averaged over " << sweeps_ << " individual spectra)." << endl;

    return spectra;
}


scalarField spectralMethodsFFTBased::frequencies
(
    const scalar& deltaT
)
{
    scalarField res(bins_/2, 0);

    forAll (res, freqi)
    {
        res[freqi] = (freqi + 1)/(deltaT*bins_);
    }

    return res;
}


Field<complex> spectralMethodsFFTBased::fft
(
    const scalarField& input
)
{
    Field<complex> res(0);

    if (sweepCount_ < sweeps_)
    {
        res.setSize( bins_/2 );

        res = computeFourierTransform( input );

        sweepCount_++;
    }

    return res;
}


List<Field<complex> > spectralMethodsFFTBased::fft
(
    const List<scalarField>& input
)
{
    List<Field<complex> > res(0);

    if (sweepCount_ < sweeps_)
    {
        res.setSize( input.size() );

        forAll (input, inputi)
        {
            Field<complex>& r( res[inputi] );
            const scalarField& in( input[inputi] );

            r.setSize( bins_/2 );

            r = computeFourierTransform( in );
        }

        sweepCount_++;
    }

    return res;
}


void spectralMethodsFFTBased::resetSweep()
{
    sweepCount_ = 0;
}


void spectralMethodsFFTBased::initSweep
(
    const scalarField& input
)
{
    resetSweep();

    sweeps_ = ( input.size() - bins_ )/step_ + 1;

    if (sweeps_ <= 0)
    {
        FatalErrorIn
        (
            "void spectralMethodsFFTBased::initSweep( const scalarField& input)"
        )
        << "The input data set is too short relative to the window size"
        << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
