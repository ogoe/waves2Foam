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
    Foam::waveTheories::bichromaticSecond

Description
    Implementation of a first order stokes theory for bichromatic surface
    water waves. For a theoretical description e.g.

        Madsen, P.A. and Fuhrman, D.R. (2006)
        Third-Order Theory for Bichromatic Bi-Directional Water Waves
        Journal of Fluid Mechanics, 557, 369-397

    Implemented without the presence of a net current and the amplitudes for
    COS and SIN are put together into one amplitude and a phase shift. This
    impacts some of the coefficients, which become identical, e.g. Anm+ = Anm-

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
    bichromaticSecond.C

Author
    Niels Gjoel Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef bichromaticSecond_H
#define bichromaticSecond_H

#include "waveTheory.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

/*---------------------------------------------------------------------------*\
                       Class bichromaticSecond Declaration
\*---------------------------------------------------------------------------*/

class bichromaticSecond
:
    public waveTheory
{

protected:

    // Protected data
        scalar Hn_;
        scalar Hm_;
        scalar h_;
        scalar omega1n_;
        scalar omega1m_;
        scalar periodn_;
        scalar periodm_;
        scalar phin_;
        scalar phim_;
        vector kn_;
        vector km_;
        scalar Tsoft_;

    // Coefficients and transfer functions
        // First order solution
        scalar an_;
        scalar am_;
        scalar Fn_;
        scalar Fm_;
        scalar kappan_;
        scalar kappam_;

        // Second order solution
        scalar Anm_;
        scalar A2n_;
        scalar A2m_;

        scalar kappanmp_;
        scalar kappanmm_;

        scalar G2n_;
        scalar G2m_;
        scalar Gnmp_;
        scalar Gnmm_;

        scalar F2n_;
        scalar F2m_;
        scalar Fnmp_;
        scalar Fnmm_;

        vector dirn_;
        vector dirm_;
        vector dirnmp_;
        vector dirnmm_;

    // Protected member functions

        void printCoeffs
        (
            const word,
            const scalar,
            const scalar
        ) const;

        void setCoeffs();

        scalar betanm
        (
            const scalar&,
            const scalar&,
            const scalar&
        ) const;

        scalar Lambda2
        (
            const scalar&,
            const vector&,
            const scalar&,
            const scalar&,
            const vector&,
            const scalar&,
            const scalar&
        ) const;

        scalar Gamma2
        (
            const scalar&,
            const vector&,
            const scalar&,
            const scalar&,
            const vector&,
            const scalar&,
            const scalar&
        ) const;

        scalar factor(const scalar&) const;

        scalar argn( const point&, const scalar&) const;

        scalar argm( const point&, const scalar&) const;

public:

    //- Runtime type information
    TypeName("bichromaticSecond");

    // Constructors

        //- from components
        bichromaticSecond
        (
            const word&,
            const fvMesh& mesh_
        );


    // Destructor

        ~bichromaticSecond()
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
