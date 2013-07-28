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
    Foam::interpolateVelocityProbes

Description
    Reads velocity probes and writes each of them into an OF-formatted
    data file.

    The data is interpolated onto an equidistant time axis.

SourceFiles
    interpolateVelocityProbes.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef interpolateVelocityProbes_H
#define interpolateVelocityProbes_H

#include "postProcessingWaves.H"
#include "rawVelocityProbes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class interpolateVelocityProbes
:
    public postProcessingWaves
{

protected:

    // Protected data

    // Protected member functions

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        interpolateVelocityProbes(const interpolateVelocityProbes&);

        //- Disallow default bitwise assignment
        void operator=(const interpolateVelocityProbes&);

    // Private member data

public:

    //- Runtime type information
    TypeName("interpolateVelocityProbes");

    // Constructors

        //- Construct from components
        interpolateVelocityProbes
        (
            const Time&,
            const dictionary&,
            const word&
        );


    // Destructor

        virtual ~interpolateVelocityProbes();

    // Member Functions
        virtual void evaluate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
