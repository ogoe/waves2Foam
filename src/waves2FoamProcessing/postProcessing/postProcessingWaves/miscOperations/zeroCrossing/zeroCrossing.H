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
    Foam::zeroCrossing

Description
    Perform an ensemble average of a given data stream.

SourceFiles
    zeroCrossing.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef zeroCrossing_H
#define zeroCrossing_H

#include "postProcessingWaves.H"

#include "spectralMethodsFFTBased.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class zeroCrossing
:
    public postProcessingWaves
{

protected:

    // Protected data
        IOdictionary dataDict_;

        labelList indices_;
        scalar deltaT_;

        word crossingType_;

    // Protected member functions

private:

    // Private Member Functions
        void evaluateScalar();

        void writeScalar
        (
            const List<scalarField>&,
            const List<scalarField>&,
            const List<scalarField>&,
            const List<scalarField>&,
            const word
        ) const;

        void crossingAnalysis
        (
            const scalarField&,
            const List<scalarField>&,
            List<scalarField>&,
            List<scalarField>&,
            List<scalarField>&,
            List<scalarField>&,
            const label
        ) const;


        //- Disallow default bitwise copy construct
        zeroCrossing(const zeroCrossing&);

        //- Disallow default bitwise assignment
        void operator=(const zeroCrossing&);

    // Private member data

public:

    //- Runtime type information
    TypeName("zeroCrossing");

    // Constructors

        //- Construct from components
        zeroCrossing
        (
            const Time&,
            const dictionary&,
            const word&
        );


    // Destructor

        virtual ~zeroCrossing();

    // Member Functions
        virtual void evaluate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
