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
    Foam::setWaveField

Description
    Implementation of a method for setting an arbitrary wave theory as initial
    conditions in the computational domain. Change in bed level is not taken
    into consideration.

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
    setWaveField.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef setWaveField_H
#define setWaveField_H

//#include "wallFvPatch.H"

#include "convexPolyhedral.H"
#include "waveTheory.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class setWaveField Declaration
\*---------------------------------------------------------------------------*/

class setWaveField
:
    public convexPolyhedral
{
private:

    // Private data

        const fvMesh& mesh_;

        volVectorField& U_;
        volScalarField& alpha_;
        volScalarField& p_;

        word name_;
        autoPtr<waveTheories::waveTheory> waveProps_;


    // Private member functions

        //- Returns a scalarField of the signed distance to an arbitrary
        //  surface
            virtual void signedPointToSurfaceDistance
            (
                const pointField&,
                scalarField&
            );

        //- Returns a scalar of the signed distance to an arbitrary surface
            virtual scalar signedPointToSurfaceDistance
            (
                const point&
            ) const;

public:

    //- Runtime type information
    TypeName("setWaveField");

    // Constructors

        //- Construct from components
        setWaveField
        (
            const fvMesh& mesh,
            volVectorField& U,
            volScalarField& alpha,
            volScalarField& p
        );

        setWaveField
        (
            const fvMesh& mesh,
            const word,
            volVectorField& U,
            volScalarField& alpha,
            volScalarField& p
        );

    // Member Functions

        void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
