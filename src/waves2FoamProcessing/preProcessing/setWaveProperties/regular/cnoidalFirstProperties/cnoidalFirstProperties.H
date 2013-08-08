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
    Foam::cnoidalFirstProperties

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
    cnoidalFirstProperties.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef cnoidalFirstProperties_H
#define cnoidalFirstProperties_H

#include "setWaveProperties.H"

namespace Foam
{

class cnoidalFirstProperties
:
    public setWaveProperties
{
private:

    scalar T_;
    scalar d_;
    scalar H_;

    double solve( );

public:

    //- Runtime type information
        TypeName("cnoidalFirstProperties");

        cnoidalFirstProperties
        (
            const Time&,
            dictionary&,
            bool
        );

    // Method

        virtual void set( Ostream&);
};

struct cnoidalFirstParams
{
    double depth_;
    double height_;
    double period_;
    double g_;
};

double lowerMBound_f ( double m, void *params);

double cnoidalFirst_f ( double m, void *params);

}

#endif

// ************************************************************************* //
