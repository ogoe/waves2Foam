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
    Foam::powerSpectraLS

Description
    This naming could be a bit misleading, however, to have the same naming
    independently from irregular and regular signal, both are called
    powerSpectra*.

    This function, however, returns the amplitudes of the cosine and sine
    contributions, which make up the signal. The output is given in the
    following format:

        1    <value> // amplitude of cos( 1*omega*t)
        1    <value> // amplitude of sin( 1*omega*t)
        2    <value> // amplitude of cos( 2*omega*t)
        2    <value> // amplitude of sin( 2*omega*t)
        ..... etc .....
        0    <value> // average value

SourceFiles
    powerSpectraLS.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef powerSpectraLS_H
#define powerSpectraLS_H

#include "postProcessingWaves.H"

#include "spectralMethodsLeastSquaresBased.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class powerSpectraLS
:
    public postProcessingWaves
{

protected:

    // Protected data
        IOdictionary dataDict_;

        labelList indices_;

        label N_;

        scalar period_;



    // Protected member functions

private:

    // Private Member Functions
        void evaluateScalar();

        void evaluateVector();

        void writeScalar
        (
            const scalarField&,
            const List<scalarField>&
        );

        void writeVector
        (
            const scalarField&,
            const List<vectorField>&
        );

        //- Disallow default bitwise copy construct
        powerSpectraLS(const powerSpectraLS&);

        //- Disallow default bitwise assignment
        void operator=(const powerSpectraLS&);

    // Private member data

public:

    //- Runtime type information
    TypeName("powerSpectraLS");

    // Constructors

        //- Construct from components
        powerSpectraLS
        (
            const Time&,
            const dictionary&,
            const word&
        );


    // Destructor

        virtual ~powerSpectraLS();

    // Member Functions
        virtual void evaluate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
