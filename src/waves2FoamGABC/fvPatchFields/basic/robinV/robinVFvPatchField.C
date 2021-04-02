/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "robinVFvPatchField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::robinVFvPatchField<Type>::robinVFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),

    sourceValue_(p.size()),
    fixedCoeffs_(p.size()),
    gradCoeffs_(p.size())
{}


template<class Type>
Foam::robinVFvPatchField<Type>::robinVFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Type& value
)
:
    fvPatchField<Type>(p, iF, value)
{}


template<class Type>
Foam::robinVFvPatchField<Type>::robinVFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict, true),

    sourceValue_("source", dict, p.size()),
    fixedCoeffs_("fixedCoeffs", dict, p.size()),
    gradCoeffs_("gradCoeffs", dict, p.size())
{}


template<class Type>
Foam::robinVFvPatchField<Type>::robinVFvPatchField
(
    const robinVFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),

    sourceValue_(ptf.sourceValue_, mapper),
    fixedCoeffs_(ptf.fixedCoeffs_, mapper),
    gradCoeffs_(ptf.gradCoeffs_, mapper)
{
    if (notNull(iF) && mapper.hasUnmapped())
    {
        WarningInFunction
            << "On field " << iF.name() << " patch " << p.name()
            << " patchField " << this->type()
            << " : mapper does not map all values." << nl
            << "    To avoid this warning fully specify the mapping in derived"
            << " patch fields." << endl;
    }
}


template<class Type>
Foam::robinVFvPatchField<Type>::robinVFvPatchField
(
    const robinVFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),

    sourceValue_(ptf.sourceValue_),
    fixedCoeffs_(ptf.fixedCoeffs_),
    gradCoeffs_(ptf.gradCoeffs_)
{}


template<class Type>
Foam::robinVFvPatchField<Type>::robinVFvPatchField
(
    const robinVFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),

    sourceValue_(ptf.sourceValue_),
    fixedCoeffs_(ptf.fixedCoeffs_),
    gradCoeffs_(ptf.gradCoeffs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::robinVFvPatchField<Type>::snGrad() const
{
    // Prepare result field
    tmp<Field<Type> > tres(new Field<Type>);

#if EXTBRANCH==1
    Field<Type>& res(tres());
#elif OFPLUSBRANCH==1
    #if OFVERSION < 1606
        Field<Type>& res(tres());
    #else
        Field<Type>& res(tres.ref());
    #endif
#else
    #if OFVERSION<400
        Field<Type>& res(tres());
    #else
        Field<Type>& res(tres.ref());
    #endif
#endif

    res.setSize(fixedCoeffs_.size());

    // Get the internal field
    Field<Type> intField = this->patchInternalField()();
    scalarField dc = this->patch().deltaCoeffs();

    // Loop over all components
    for (int cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        const scalarField a  = fixedCoeffs_.component(cmpt);
        const scalarField b  = gradCoeffs_.component(cmpt);
        const scalarField r  = sourceValue_.component(cmpt);
        const scalarField IF = intField.component(cmpt);

        res.replace(cmpt, dc*(r - a*IF)/stabilise(a + b*dc, SMALL));
    }

    return tres;
}


template<class Type>
void Foam::robinVFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // Prepare result field
    Field<Type> res(fixedCoeffs_.size());

    // Get the internal field
    Field<Type> intField = this->patchInternalField()();
    scalarField dc = this->patch().deltaCoeffs();

    // Loop over all components
    for (int cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        const scalarField a  = fixedCoeffs_.component(cmpt);
        const scalarField b  = gradCoeffs_.component(cmpt);
        const scalarField r  = sourceValue_.component(cmpt);
        const scalarField IF = intField.component(cmpt);

        res.replace(cmpt, (b*dc*IF + r)/stabilise(a + b*dc, SMALL));
    }

    Field<Type>::operator=(res);

    fvPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::robinVFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    // Prepare result field
    tmp<Field<Type> > tres(new Field<Type>);

#if EXTBRANCH==1
    Field<Type>& res(tres());
#elif OFPLUSBRANCH==1
    #if OFVERSION < 1606
        Field<Type>& res(tres());
    #else
        Field<Type>& res(tres.ref());
    #endif
#else
    #if OFVERSION<400
        Field<Type>& res(tres());
    #else
        Field<Type>& res(tres.ref());
    #endif
#endif

    res.setSize(fixedCoeffs_.size());

    // Loop over all components
    scalarField dc = this->patch().deltaCoeffs();

    for (int cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        const scalarField a  = fixedCoeffs_.component(cmpt);
        const scalarField b  = gradCoeffs_.component(cmpt);

        res.replace(cmpt, b*dc/stabilise(a + b*dc, SMALL));
    }

    return tres;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::robinVFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    // Prepare result field
    tmp<Field<Type> > tres(new Field<Type>);

#if EXTBRANCH==1
    Field<Type>& res(tres());
#elif OFPLUSBRANCH==1
    #if OFVERSION < 1606
        Field<Type>& res(tres());
    #else
        Field<Type>& res(tres.ref());
    #endif
#else
    #if OFVERSION<400
        Field<Type>& res(tres());
    #else
        Field<Type>& res(tres.ref());
    #endif
#endif

    res.setSize(fixedCoeffs_.size());

    // Loop over all components
    scalarField dc = this->patch().deltaCoeffs();

    for (int cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        const scalarField a  = fixedCoeffs_.component(cmpt);
        const scalarField b  = gradCoeffs_.component(cmpt);
        const scalarField r  = sourceValue_.component(cmpt);

        res.replace(cmpt, r/stabilise(a + b*dc, SMALL));
    }

    return tres;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::robinVFvPatchField<Type>::gradientInternalCoeffs() const
{
    // Prepare result field
    tmp<Field<Type> > tres(new Field<Type>);

#if EXTBRANCH==1
    Field<Type>& res(tres());
#elif OFPLUSBRANCH==1
    #if OFVERSION < 1606
        Field<Type>& res(tres());
    #else
        Field<Type>& res(tres.ref());
    #endif
#else
    #if OFVERSION<400
        Field<Type>& res(tres());
    #else
        Field<Type>& res(tres.ref());
    #endif
#endif

    res.setSize(fixedCoeffs_.size());

    // Loop over all components
    scalarField dc = this->patch().deltaCoeffs();

    for (int cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        const scalarField a  = fixedCoeffs_.component(cmpt);
        const scalarField b  = gradCoeffs_.component(cmpt);

        res.replace(cmpt, -dc*a/stabilise(a + b*dc, SMALL));
    }

    return tres;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::robinVFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    // Prepare result field
    tmp<Field<Type> > tres(new Field<Type>);

#if EXTBRANCH==1
    Field<Type>& res(tres());
    #elif OFPLUSBRANCH==1
    #if OFVERSION < 1606
        Field<Type>& res(tres());
    #else
        Field<Type>& res(tres.ref());
    #endif
#else
    #if OFVERSION<400
        Field<Type>& res(tres());
    #else
       Field<Type>& res(tres.ref());
    #endif
#endif

    res.setSize(fixedCoeffs_.size());

    // Get the internal field
    Field<Type> intField = this->patchInternalField()();

    // Loop over all components
    scalarField dc = this->patch().deltaCoeffs();

    for (int cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        const scalarField a  = fixedCoeffs_.component(cmpt);
        const scalarField b  = gradCoeffs_.component(cmpt);
        const scalarField r  = sourceValue_.component(cmpt);

        res.replace(cmpt, r*dc/stabilise(a + b*dc, SMALL));
    }

    return tres;
}


template<class Type>
void Foam::robinVFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->writeEntry("value", os);

    sourceValue_.writeEntry("source", os);
    fixedCoeffs_.writeEntry("fixedCoeffs", os);
    gradCoeffs_.writeEntry("gradCoeffs", os);
}


// ************************************************************************* //
