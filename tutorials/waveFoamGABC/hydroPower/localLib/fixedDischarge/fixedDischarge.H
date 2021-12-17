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
    Foam::waveTheories::fixedDischargeTrial

Description

SourceFiles
    fixedDischargeTrial.C

Author
    Niels Gjoel Jacobsen.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef fixedDischargeTrial_H
#define fixedDischargeTrial_H

#include "waveTheory.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

/*---------------------------------------------------------------------------*\
                           Class fixedDischargeTrial Declaration
\*---------------------------------------------------------------------------*/

class fixedDischargeTrial
:
    public waveTheory
{

protected:

    // Protected data
        const fvMesh& fvMesh_;
        word boundaryName_;

        scalar Tsoft_;
        scalar Qm3s_;

        scalar bottomLevel_;
        scalar width_;

        mutable label timeIndex_;
        mutable scalar localSeaLevel_;

        scalar filterTime_;
        scalar filterDt_;
        mutable scalar filterLastSample_;

        mutable scalarField surfaceTime_;
        mutable scalarField surfaceElevation_;

    // Protected member functions
        void updateLocalSeaLevel() const;

        scalar factor(const scalar&) const;
public:

    //- Runtime type information
    TypeName("fixedDischargeTrial");

    // Constructors

        //- from components
        fixedDischargeTrial
        (
            const word&,
            const fvMesh& mesh_
        );


    // Destructor

        ~fixedDischargeTrial()
        {}


    // Member Functions
        virtual Switch conflictTSoftInitialise() const
        {
        	return true;
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
