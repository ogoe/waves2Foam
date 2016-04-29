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
    Foam::probeGauges

Description
    Writing in input file for probes based on a given point distribution.

SourceFiles
    probeGauges.C

Author
    Niels Gjøl Jacobsen, Deltares.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef probeGauges_H
#define probeGauges_H

#include "fvMesh.H"
#include "fvCFD.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "OFstream.H"
#include "uniformDimensionedFields.H"

#include "pointDistributions.H"
#include "crossVersionCompatibility.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class porosityCoefficient Declaration
\*---------------------------------------------------------------------------*/

class probeGauges
{
private:

    // Private data
        const fvMesh& mesh_;
        const dictionary& gaugeDict_;

    // Private Member Functions
        void writeVTKFormat
        (
            const word&,
            const pointField&
        );

        //- Disallow default bitwise copy construct
        probeGauges(const probeGauges&);

        //- Disallow default bitwise assignment
        void operator=(const probeGauges&);

public:

    //- Runtime type information
    TypeName("probeGauges");

    // Constructors

        //- Construct from components
        probeGauges
        (
            const fvMesh&,
            const dictionary&
        );

    // Destructor

        virtual ~probeGauges() {};

    // Member Functions

        // Access
        void evaluate(const word&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
