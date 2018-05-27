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
    Foam::waveTheories::externalSource

Description


SourceFiles
    externalSource.C

Author
    Niels Gjoel Jacobsen, Deltares.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef externalSource_H
#define externalSource_H

#include "waveTheory.H"
#include "externalWaveForcing.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

/*---------------------------------------------------------------------------*\
                           Class externalSource Declaration
\*---------------------------------------------------------------------------*/

class externalSource
:
    public waveTheory
{

protected:

    // Protected data
        const externalWaveForcing& external_;

    // Protected member functions

        scalar factor(const scalar&) const;
public:

    //- Runtime type information
    TypeName("externalSource");

    // Constructors

        //- from components
        externalSource
        (
            const word&,
            const fvMesh& mesh_
        );


    // Destructor

        ~externalSource()
        {}


    // Member Functions
        virtual Switch conflictTSoftInitialise() const
        {
            return true;
        };

        void printCoeffs();

        virtual scalar eta
        (
            const point&,
            const scalar&
        ) const;

//        virtual scalar ddxPd
//        (
//            const point&,
//            const scalar&,
//            const vector&
//        ) const;

        virtual scalar pExcess
        (
            const point&,
            const scalar&
        ) const;

        bool implementPressure() const
        {
            return external_.implementPressure();
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
