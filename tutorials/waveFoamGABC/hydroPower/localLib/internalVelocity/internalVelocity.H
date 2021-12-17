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
    Foam::waveTheories::internalVelocity

Description


SourceFiles
    internalVelocity.C

Author
    Niels Gjoel Jacobsen.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef internalVelocity_H
#define internalVelocity_H

#include "waveTheory.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

/*---------------------------------------------------------------------------*\
                           Class internalVelocity Declaration
\*---------------------------------------------------------------------------*/

class internalVelocity
:
    public waveTheory
{

protected:

    // Protected data
        const fvMesh& fvMesh_;
        word boundaryName_;

        scalar localSeaLevel_;

        scalar Tsoft_;

        mutable label timeIndex_;
        mutable scalarField phiBoundary_;

    // Protected member functions
        const scalarField& boundaryFlux() const;

        scalar factor(const scalar&) const;
public:

    //- Runtime type information
    TypeName("internalVelocityTrial");

    // Constructors

        //- from components
        internalVelocity
        (
            const word&,
            const fvMesh& mesh_
        );


    // Destructor

        ~internalVelocity()
        {}


    // Member Functions
        virtual Switch conflictTSoftInitialise() const
        {
        	return false;
        };

        virtual bool isInitialiser() const
        {
            return false;
        };

        void printCoeffs();

        virtual scalar eta
        (
            const point&,
            const scalar&
        ) const;

        virtual scalar pExcess
        (
            const point&,
            const scalar&
        ) const;

        bool implementPressure() const
        {
            return true;
        };

        virtual vector U
        (
            const point&,
            const scalar&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
