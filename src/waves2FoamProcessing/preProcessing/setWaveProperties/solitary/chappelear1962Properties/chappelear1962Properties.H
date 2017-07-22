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
    Foam::chappelear1962Properties

Description

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
    chappelear1962Properties.C

Author
    Niels G. Jacobsen, Deltares.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef chappelear1962Properties_H
#define chappelear1962Properties_H

#include "setWaveProperties.H"

namespace Foam
{

class chappelear1962Properties
:
    public setWaveProperties
{
private:

	scalar G1(const scalar, const scalar) const;

	scalar G2(const scalar, const scalar, const scalar) const;

	scalarField Jacobi(const scalar, const scalar) const;

public:

        //- Runtime type information
        TypeName("chappelear1962Properties");

        chappelear1962Properties
        (
            const Time&,
            dictionary&,
            bool
        );

        // Method

        virtual void set(Ostream&);
};

}

#endif

// ************************************************************************* //
