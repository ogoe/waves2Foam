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
    Foam::waveTheories::chappelear1962

Description
    Implementation of the solitary wave theory described in Chappelear (1962).
    Shallow-Water Waves, Journal of Geophysical Research.

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
    chappelear1962.C

Author
    Niels G. Jacobsen, Deltares.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef chappelear1962_H
#define chappelear1962_H

#include "waveTheory.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

/*---------------------------------------------------------------------------*\
                         Class chappelear1962 Declaration
\*---------------------------------------------------------------------------*/

class chappelear1962
:
    public waveTheory
{

protected:

    // Protected data
        // Read from waveProperties
        scalar H_;
        scalar h_;
        vector propagationDirection_;
        point  x0_;
        scalar L1_;
        scalar L3_;

        // Bernoullis constant - estimate
        scalar R_;

        // Computed
        scalar dL_;
        scalar kappa_;
        scalar G_;
        scalar c_;

    // Protected member functions

        scalar shapeT(const point&, const scalar&) const;

        scalar shapeS(const point&, const scalar&) const;

        scalar factor(const scalar& ) const;
public:

    //- Runtime type information
    TypeName("chappelear1962");

    // Constructors

        //- from components
        chappelear1962
        (
            const word&,
            const fvMesh& mesh_
        );


    // Destructor

        ~chappelear1962()
        {}


    // Member Functions
        virtual Switch conflictTSoftInitialise() const
        {
            return false;
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

        virtual scalar pExcess
        (
            const point&,
            const scalar&
        ) const;

        bool implementPressure() const
        {
            return false;
        };

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
