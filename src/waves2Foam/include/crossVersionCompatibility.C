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

\*---------------------------------------------------------------------------*/

#include "crossVersionCompatibility.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace waves2Foam
{

word pName()
{
#if EXTBRANCH==1
    return "pd";
#elif OFPLUSBRANCH==1
    return "p_rgh";
#else
    #if OFVERSION<170
        return "pd";
    #else
        return "p_rgh";
    #endif
#endif
}


word aName()
{
#if EXTBRANCH==1
    return "alpha1";
#elif OFPLUSBRANCH==1
    return "alpha.water";
#else
    #if OFVERSION<230
        return "alpha1";
    #else
        return "alpha.water";
    #endif
#endif
}


word waterPhase()
{
#if EXTBRANCH==1
    return "phase1";
#elif OFPLUSBRANCH==1
    return "water";
#else
    #if OFVERSION<230
        return "phase1";
    #else
        return "water";
    #endif
#endif
}


word airPhase()
{
#if EXTBRANCH==1
    return "phase2";
#elif OFPLUSBRANCH==1
    return "air";
#else
    #if OFVERSION<230
        return "phase2";
    #else
        return "air";
    #endif
#endif

}


word rAUfName()
{
#if EXTBRANCH==1
	return "Dp";
#elif OFPLUSBRANCH==1
	return "rAUf";
#elif 290<OFVERSION
	return "rAUf";
#else
	notImplemented("word rAUfName() is not implemented for this OpenFoam version.");
#endif
}


} // End namespace waves2Foam

} // End namespace Foam

// ************************************************************************* //
