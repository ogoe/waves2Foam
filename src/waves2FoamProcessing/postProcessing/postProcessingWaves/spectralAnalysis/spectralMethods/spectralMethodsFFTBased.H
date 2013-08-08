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

Class
    Foam::spectralMethodsFFTBased

Description
    Helper functions for FFT-based analysis.

SourceFiles
    spectralMethodsFFTBased.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef spectralMethodsFFTBased_H
#define spectralMethodsFFTBased_H

#include "fvCFD.H"
#include "complexExp.H"

#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class spectralMethodsFFTBased
{

private:
    // Private member data
        const dictionary& dict_;

        label bins_;
        label step_;
        label sweeps_;
        label sweepCount_;

        // Initialise and allocate GSL-FFT
        gsl_fft_real_wavetable *real;
        gsl_fft_real_workspace *work;

    // Private Member Functions
        void checkBins();

        void powerSpectrum
        (
            const scalarField&,
            const scalar&,
            scalarField&
        );

        void powerSpectrum
        (
            const vectorField&,
            const scalar&,
            vectorField&
        );

        Field<complex> computeFourierTransform
        (
            const scalarField&
        );

        //- Disallow default bitwise copy construct
        spectralMethodsFFTBased(const spectralMethodsFFTBased&);

        //- Disallow default bitwise assignment
        void operator=(const spectralMethodsFFTBased&);

public:

    //- Runtime type information
    TypeName("spectralMethodsFFTBased");

    // Constructors

        //- Construct from components
        spectralMethodsFFTBased
        (
            const Time&,
            const dictionary&
        );

    // Destructor

        virtual ~spectralMethodsFFTBased();

    // Member Functions
        Field<complex> fft
        (
            const scalarField&
        );

        List<Field<complex> > fft
        (
            const List<scalarField>&
        );

        List<scalarField> powerSpectra
        (
            const List<scalarField>&,
            const scalar&
        );

        List<vectorField> powerSpectra
        (
            const List<vectorField>&,
            const scalar&
        );

        scalarField frequencies
        (
            const scalar& deltaT
        );

        void initSweep
        (
            const scalarField&
        );

        void resetSweep();

        inline label nSweeps() const { return sweeps_; };

        inline label nBins() const { return bins_; };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
