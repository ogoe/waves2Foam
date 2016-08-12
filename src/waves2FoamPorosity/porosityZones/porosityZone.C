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

#include "porosityZone.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "geometricOneField.H"

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

// adjust negative resistance values to be multiplier of max value
void Foam::porosityZone::checkNegativeResistance(dimensionedVector& resist)
{
    scalar minCmpt = min(0, cmptMin(resist.value()));

    if (minCmpt < 0)
    {
        FatalErrorIn
        (
            "Foam::porosityZone::porosityZone::adjustNegativeResistance"
            "(dimensionedVector&)"
        )   << "negative resistances! " << resist
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityZone::porosityZone
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    mesh_(mesh),
    dict_(dict),
    cellZoneID_(mesh_.cellZones().findZoneID(name)),

#if EXTBRANCH==1
    coordSys_(dict, mesh),
#elif OFPLUSBRANCH==1
    coordSys_(mesh, dict),
#else
    #if OFVERSION<230
        coordSys_(dict, mesh),
    #else
        coordSys_(mesh, dict),
    #endif
#endif

//#if OFVERSION<230 || EXTBRANCH == 1
//    coordSys_(dict, mesh),
//#else
//    coordSys_(mesh, dict),
//#endif

    porosity_( readScalar( dict_.lookup("porosity") ) ),
    addedMassCoeff_( readScalar( dict_.lookup("gammaAddedMass") ) ),
    D_("D", dimensionSet(0, -2, 0, 0, 0), tensor::zero),
    F_("F", dimensionSet(0, -1, 0, 0, 0), tensor::zero)
{
    Info<< "Creating porous zone: " << name_ << endl;

    autoPtr<Foam::porosityCoefficient> pcPtr( Foam::porosityCoefficient::New( dict ) );

    bool foundZone = (cellZoneID_ != -1);
    reduce(foundZone, orOp<bool>());

    if (!foundZone && Pstream::master())
    {
        FatalErrorIn
        (
            "Foam::porosityZone::porosityZone"
            "(const fvMesh&, const word&, const dictionary&)"
        )   << "cannot find porous cellZone " << name_
            << exit(FatalError);
    }

    // porosity
    if (porosity_ <= 0.0 || porosity_ > 1.0)
    {
        FatalIOErrorIn
        (
                "Foam::porosityZone::porosityZone"
                "(const fvMesh&, const word&, const dictionary&)",
                dict_
        )
        << "out-of-range porosity value " << porosity_
        << exit(FatalIOError);
    }

    // local-to-global transformation tensor
#if EXTBRANCH==1
    const tensor& E = coordSys_.R();
#elif OFPLUSBRANCH==1
    const tensor E = coordSys_.R().R();
#else
    #if OFVERSION<230
        const tensor& E = coordSys_.R();
    #else
        const tensor E = coordSys_.R().R();
    #endif
#endif

    dimensionedVector d( pcPtr->linearCoefficient() );

    if (D_.dimensions() != d.dimensions())
    {
        FatalIOErrorIn
        (
            "Foam::porosityZone::porosityZone"
            "(const fvMesh&, const word&, const dictionary&)",
            dict_
        )   << "incorrect dimensions for d: " << d.dimensions()
        << " should be " << D_.dimensions()
        << exit(FatalIOError);
    }

    checkNegativeResistance(d);

    D_.value().xx() = d.value().x();
    D_.value().yy() = d.value().y();
    D_.value().zz() = d.value().z();
    D_.value() = (E & D_ & E.T()).value();


    dimensionedVector f( pcPtr->quadraticCoefficient() );

    if (F_.dimensions() != f.dimensions())
    {
        FatalIOErrorIn
        (
            "Foam::porosityZone::porosityZone"
            "(const fvMesh&, const word&, const dictionary&)",
            dict_
        )   << "incorrect dimensions for f: " << f.dimensions()
        << " should be " << F_.dimensions()
        << exit(FatalIOError);
    }

    checkNegativeResistance(f);

    F_.value().xx() = 0.5 * f.value().x();
    F_.value().yy() = 0.5 * f.value().y();
    F_.value().zz() = 0.5 * f.value().z();
    F_.value() = (E & F_ & E.T()).value();

    // it is an error not to define anything
    if
    (
        magSqr(D_.value()) <= VSMALL
     && magSqr(F_.value()) <= VSMALL
    )
    {
        FatalIOErrorIn
        (
            "Foam::porosityZone::porosityZone"
            "(const fvMesh&, const word&, const dictionary&)",
            dict_
        )   << "neither powerLaw (C0/C1) "
               "nor Darcy-Forchheimer law (d/f) specified"
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porosityZone::porosity
(
    volScalarField & poro
) const
{
    if ( cellZoneID_ == -1 )
    {
        return;
    }

    const labelList& cells = mesh_.cellZones()[cellZoneID_];

    forAll( cells, celli )
    {
        poro[cells[celli]] = porosity_;
    }
}


void Foam::porosityZone::addResistance(fvVectorMatrix& UEqn) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    bool compressible = false;
    if (UEqn.dimensions() == dimensionSet(1, 1, -2, 0, 0))
    {
        compressible = true;
    }

    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();
    const vectorField& U = UEqn.psi();

    const tensor& D = D_.value();
    const tensor& F = F_.value();

    if (magSqr(D) > VSMALL || magSqr(F) > VSMALL)
    {
        if (compressible)
        {
            addViscousInertialResistance
            (
                Udiag,
                Usource,
                cells,
                V,
                mesh_.lookupObject<volScalarField>("rho"),
                mesh_.lookupObject<volScalarField>("mu"),
                U
            );
        }
        else
        {
            addViscousInertialResistance
            (
                Udiag,
                Usource,
                cells,
                V,
                geometricOneField(),
                mesh_.lookupObject<volScalarField>("nu"),
                U
            );
        }
    }
}


void Foam::porosityZone::addResistance
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU,
    bool correctAUprocBC
) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    bool compressible = false;
    if (UEqn.dimensions() == dimensionSet(1, 1, -2, 0, 0))
    {
        compressible = true;
    }

    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    const vectorField& U = UEqn.psi();

    const tensor& D = D_.value();
    const tensor& F = F_.value();

    if (magSqr(D) > VSMALL || magSqr(F) > VSMALL)
    {
        if (compressible)
        {
            addViscousInertialResistance
            (
                AU,
                cells,
                mesh_.lookupObject<volScalarField>("rho"),
                mesh_.lookupObject<volScalarField>("mu"),
                U
            );
        }
        else
        {
            addViscousInertialResistance
            (
                AU,
                cells,
                geometricOneField(),
                mesh_.lookupObject<volScalarField>("nu"),
                U
            );
        }
    }

    if (correctAUprocBC)
    {
        // Correct the boundary conditions of the tensorial diagonal to ensure
        // processor boundaries are correctly handled when AU^-1 is interpolated
        // for the pressure equation.
        AU.correctBoundaryConditions();
    }
}


void Foam::porosityZone::writeDict(Ostream& os, bool subDict) const
{
    if (subDict)
    {
        os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
        os.writeKeyword("name") << zoneName() << token::END_STATEMENT << nl;
    }
    else
    {
        os  << indent << zoneName() << nl
            << indent << token::BEGIN_BLOCK << incrIndent << nl;
    }

    if (dict_.found("note"))
    {
        os.writeKeyword("note") << string(dict_.lookup("note"))
            << token::END_STATEMENT << nl;
    }

    coordSys_.writeDict(os, true);

    if (dict_.found("porosity"))
    {
        os.writeKeyword("porosity") << porosity() << token::END_STATEMENT << nl;
    }

    // powerLaw coefficients
    if (const dictionary* dictPtr = dict_.subDictPtr("powerLaw"))
    {
        os << indent << "powerLaw";
        dictPtr->write(os);
    }

    // Darcy-Forchheimer coefficients
    if (const dictionary* dictPtr = dict_.subDictPtr("Darcy"))
    {
        os << indent << "Darcy";
        dictPtr->write(os);
    }

    os << decrIndent << indent << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const porosityZone& pZone)
{
    pZone.writeDict(os);
    return os;
}

// ************************************************************************* //
