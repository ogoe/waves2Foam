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
    Foam::interpolateSurfaceElevation

Description
    Reads probed alpha1-values and writes data for each probe in an OF-
    formatted manner.

    The data is interpolated onto an equidistant time axis.

SourceFiles
    interpolateAlphaProbes.C

Author
    Bjarne Jensen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef interpolateAlphaProbes_H
#define interpolateAlphaProbes_H

#include "postProcessingWaves.H"
#include "rawAlphaProbes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class interpolateAlphaProbes
:
    public postProcessingWaves
{

protected:

    // Protected data

    // Protected member functions

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        interpolateAlphaProbes(const interpolateAlphaProbes&);

        //- Disallow default bitwise assignment
        void operator=(const interpolateAlphaProbes&);

    // Private member data

public:

    //- Runtime type information
    TypeName("interpolateAlphaProbes");

    // Constructors

        //- Construct from components
        interpolateAlphaProbes
        (
            const Time&,
            const dictionary&,
            const word&
        );


    // Destructor

        virtual ~interpolateAlphaProbes();

    // Member Functions
        virtual void evaluate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
