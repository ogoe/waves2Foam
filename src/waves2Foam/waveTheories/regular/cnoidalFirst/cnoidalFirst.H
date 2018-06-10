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
    Foam::waveTheories::cnoidalFirst

Description
    Implementation of a first order cnoidal theory for monochromatic surface
    water waves. Implemented following the description by

        Hydrodynamics of Coastal Regions
        Ib A. Svendsen and Ivar G. Jonsson
        1982
        Den private ingenioerfond
        Technical University of Denmark

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
    cnoidalFirst.C

Author
    Niels Gjoel Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef cnoidalFirst_H
#define cnoidalFirst_H

#include "waveTheory.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

/*---------------------------------------------------------------------------*\
                           Class cnoidalFirst Declaration
\*---------------------------------------------------------------------------*/

class cnoidalFirst
:
    public waveTheory
{

protected:

    // Protected data
        scalar H_;
        scalar h_;
        scalar omega_;
        vector propagationDirection_;
        scalar m_;
        scalar length_;
        scalar celerity_;
        scalar Tsoft_;

        scalar period_;
        scalar Kelliptic_;
        scalar etaMin_;
    // Protected member functions

        scalar factor(const scalar&) const;

    // Helper functions
        scalar argument
        (
            const point&,
            const scalar&
        ) const;

        scalar eta_x
        (
            const scalar&,
            const scalar&,
            const scalar&
        ) const;

        scalar eta_xx
        (
            const scalar&,
            const scalar&,
            const scalar&
        ) const;

        scalar eta_xxx
        (
            const scalar&,
            const scalar&,
            const scalar&
        ) const;

public:

    //- Runtime type information
    TypeName("cnoidalFirst");

    // Constructors

        //- from components
        cnoidalFirst
        (
            const word&,
            const fvMesh& mesh_
        );


    // Destructor

        ~cnoidalFirst()
        {}


    // Member Functions
        virtual Switch conflictTSoftInitialise() const
        {
            if (Tsoft_ > 0)
            {
                return true;
            }
            else
            {
                return false;
            }
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
