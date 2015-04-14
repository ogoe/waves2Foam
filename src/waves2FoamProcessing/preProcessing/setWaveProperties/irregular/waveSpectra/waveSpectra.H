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
    Foam::waveSpectra

Description
    An abstract base class for a pre-processing tool for setting the wave
    properties in the dictionary waveProperties.

    A description of the general wave framework is given in

    @article { jacobsenFuhrmanFredsoe2011,
        Author = {Jacobsen, N G and Fuhrman, D R and Freds\o{}e, J},
        title = {{A Wave Generation Toolbox for the Open-Source CFD Library: OpenFoam\textregistered{}}},
        Journal = {{Int. J. for Numer. Meth. Fluids}},
        Year = {2012},
        Volume = {70},
        Number = {9},
        Pages = {1073-1088},
        DOI = {{10.1002/fld.2726}},
    }

SourceFiles
    waveSpectra.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef waveSpectra_H
#define waveSpectra_H

#include "fvMesh.H"
#include "volFields.H"
#include "IOdictionary.H"

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "uniformDimensionedFields.H"

#include "setWaveProperties.H"
#include "equidistantFrequencyAxis.H"
#include "phases.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class waveSpectra Declaration
\*---------------------------------------------------------------------------*/

class waveSpectra
{

protected:

    // Protected data
        const Time& rT_;

        dictionary& dict_;

        scalarField& amp_;
        scalarField& freq_;
        scalarField& phi_;
        vectorField& k_;

        scalar G_;

        scalar PI_;

        autoPtr<Foam::phases> phases_;


    // Protected member functions
        void writeSpectrum
        (
            Ostream&,
            const scalarField&,
            const scalarField&
        ) const;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        waveSpectra(const waveSpectra&);

        //- Disallow default bitwise assignment
        void operator=(const waveSpectra&);

public:

    //- Runtime type information
    TypeName("waveSpectra");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            waveSpectra,
            waveSpectra,
            (
                const Time& rT,
                dictionary& dict,
                scalarField& amp,
                scalarField& freq,
                scalarField& phi,
                vectorField& k
            ),
            (rT, dict, amp, freq, phi, k)
        );

    // Constructors

        //- Construct from components
        waveSpectra
        (
            const Time&,
            dictionary&,
            scalarField&,
            scalarField&,
            scalarField&,
            vectorField&
        );

    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<waveSpectra> New
        (
            const Time&,
            dictionary&,
            scalarField&,
            scalarField&,
            scalarField&,
            vectorField&
        );

    // Destructor

        virtual ~waveSpectra();

    // Member functions

        virtual void set(Ostream&) = 0;

        virtual wordList list() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
