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
    Foam::concatenatedDistribution

Description
    A simple concatenated distribution with stretching.

SourceFiles
    concatenatedDistribution.C

Author
    Niels Gjoel Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef concatenatedDistribution_H
#define concatenatedDistribution_H

#include "pointDistributions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class concatenatedDistribution
:
    public pointDistributions
{

protected:

    // Protected data
        pointField ppE_;
        pointField ppS_;

    // Protected member functions

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        concatenatedDistribution(const concatenatedDistribution&);

        //- Disallow default bitwise assignment
        void operator=(const concatenatedDistribution&);

    // Private member data

public:

    //- Runtime type information
    TypeName("concatenatedDistribution");

    // Constructors

        //- Construct from components
        concatenatedDistribution
        (
//            const fvMesh&,
            const dictionary&
        );


    // Destructor

        virtual ~concatenatedDistribution();

    // Member Functions
        virtual pointField evaluateStart();

        virtual pointField evaluateEnd();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
