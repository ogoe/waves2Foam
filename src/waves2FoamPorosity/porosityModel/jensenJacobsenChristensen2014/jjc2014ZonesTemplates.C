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

\*----------------------------------------------------------------------------*/

#include "jjc2014Zones.H"
#include "volFields.H"
#include "fvMatrix.H"
#include "fvm.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::jjc2014Zones::modifyDdt(fvMatrix<Type>& m) const
{
    forAll(*this, i)
    {
        operator[](i).modifyDdt(m);
    }
}


// * * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> >
Foam::jjc2014Zones::ddt
(
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tres = fvm::ddt(vf);
#if EXTBRANCH==1
    modifyDdt(tres());
#elif OFPLUSBRANCH==1
    #if OFVERSION < 1606
        modifyDdt(tres());
    #else
        modifyDdt(tres.ref());
    #endif
#else
    #if OFVERSION<400
        modifyDdt(tres());
    #else
        modifyDdt(tres.ref());
    #endif
#endif
    
    return tres;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> >
Foam::jjc2014Zones::ddt
(
    const geometricOneField&,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tres = fvm::ddt(vf);
#if EXTBRANCH==1
    modifyDdt(tres());
#elif OFPLUSBRANCH==1
    #if OFVERSION < 1606
        modifyDdt(tres());
    #else
        modifyDdt(tres.ref());
    #endif
#else
    #if OFVERSION<400
        modifyDdt(tres());
    #else
        modifyDdt(tres.ref());
    #endif
#endif

    return tres;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> >
Foam::jjc2014Zones::ddt
(
    const dimensionedScalar& rho,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tres = fvm::ddt(rho,vf);
#if EXTBRANCH==1
    modifyDdt(tres());
#elif OFPLUSBRANCH==1
    #if OFVERSION < 1606
        modifyDdt(tres());
    #else
        modifyDdt(tres.ref());
    #endif
#else
    #if OFVERSION<400
        modifyDdt(tres());
    #else
        modifyDdt(tres.ref());
    #endif
#endif

    return tres;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> >
Foam::jjc2014Zones::ddt
(
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tres = fvm::ddt(rho,vf);
#if EXTBRANCH==1
    modifyDdt(tres());
#elif OFPLUSBRANCH==1
    #if OFVERSION < 1606
        modifyDdt(tres());
    #else
        modifyDdt(tres.ref());
    #endif
#else
    #if OFVERSION<400
        modifyDdt(tres());
    #else
        modifyDdt(tres.ref());
    #endif
#endif

    return tres;
}

// ************************************************************************* //
