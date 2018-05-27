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
    Foam::waveTheories::combinedWaves

Description
    Allow for a arbitrary combination of the other wave theories. Not tested,
    but should also be able to work recursively.

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
    combinedWaves.C

\*---------------------------------------------------------------------------*/

#ifndef combinedWaves_H
#define combinedWaves_H

#include "waveTheory.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

/*---------------------------------------------------------------------------*\
                           Class combinedWaves Declaration
\*---------------------------------------------------------------------------*/

class combinedWaves
:
    public waveTheory
{

protected:

    // Protected data
    wordList combinedWavesNames_;

    List<autoPtr<waveTheories::waveTheory> > combinedWavesPtr_;
    // Protected member functions

        scalar factor(const scalar&) const;
public:

    //- Runtime type information
    TypeName("combinedWaves");

    // Constructors

        //- from components
        combinedWaves
        (
            const word&,
            const fvMesh& mesh_
        );


    // Destructor

        ~combinedWaves()
        {}


    // Member Functions
        virtual Switch conflictTSoftInitialise() const
        {
            Switch res = false;

            forAll (combinedWavesPtr_, wavei)
            {
                res = combinedWavesPtr_[wavei]->conflictTSoftInitialise();

                if (res)
                {
                    break;
                }
            }

            return res;
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
            forAll (combinedWavesPtr_, wavei)
            {
                if (~combinedWavesPtr_[wavei]->implementPressure())
                {
                    return false;
                }
            }
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
