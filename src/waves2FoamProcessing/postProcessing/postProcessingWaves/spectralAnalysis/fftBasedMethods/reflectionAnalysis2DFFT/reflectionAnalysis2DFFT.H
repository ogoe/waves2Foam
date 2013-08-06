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
    Foam::reflectionAnalysis2DFFT

Description
    Compute the incident and reflected spectra for a given set of surface
    elevation gauges. The method is to be used on horizontal beds and it
    uses the linear dispersion relation to estimate the wave numbers.

    The method for decomposing the spectra is that of Zelt and Skjelbreia.

SourceFiles
    reflectionAnalysis2DFFT.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef reflectionAnalysis2DFFT_H
#define reflectionAnalysis2DFFT_H

#include "postProcessingWaves.H"

#include "stokesFirstProperties.H"
#include "spectralMethodsFFTBased.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class reflectionAnalysis2DFFT
:
    public postProcessingWaves
{

protected:

    // Protected data
        IOdictionary dataDict_;

        labelList indices_;
        scalar deltaT_;
        scalar depth_;

        word coordName_;
        label writeIndex_;

        scalarField X_;

    // Protected member functions

private:

    // Private Member Functions
        scalarField linearWaveNumbers
        (
            const scalarField&
        );

        void writeReflectionIncident
        (
            const scalarField&,
            const scalarField&,
            const scalarField&,
            const scalarField&
        );

        void decomposeAmplitudes
        (
            const scalarField&,
            const List<Field<complex> >&,
            Field<complex>&,
            Field<complex>&,
            scalarField&
        );

        //- Disallow default bitwise copy construct
        reflectionAnalysis2DFFT(const reflectionAnalysis2DFFT&);

        //- Disallow default bitwise assignment
        void operator=(const reflectionAnalysis2DFFT&);

    // Private member data

public:

    //- Runtime type information
    TypeName("reflectionAnalysis2DFFT");

    // Constructors

        //- Construct from components
        reflectionAnalysis2DFFT
        (
            const Time&,
            const dictionary&,
            const word&
        );

    // Destructor

        virtual ~reflectionAnalysis2DFFT();

    // Member Functions
        virtual void evaluate();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
