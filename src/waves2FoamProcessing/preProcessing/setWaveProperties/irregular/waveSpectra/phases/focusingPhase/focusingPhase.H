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
    Foam::focusingPhase

Description
    Set phases such that the irregular wave train focuses at a given time at
    a given location.

SourceFiles
    focusingPhase.C

Author
    Niels Gjøl Jacobsen, Deltares.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef focusingPhase_H
#define focusingPhase_H

#include "phases.H"
#include <time.h>

namespace Foam
{

class focusingPhase
:
    public phases
{
private:

        scalar focusTime_;
        vector focusPoint_;

public:

        //- Runtime type information
        TypeName("focusingPhase");

        focusingPhase
        (
            const Time&,
            dictionary&
        );

        virtual scalar phase(const scalar&, const vector&);
};

}

#endif

// ************************************************************************* //
