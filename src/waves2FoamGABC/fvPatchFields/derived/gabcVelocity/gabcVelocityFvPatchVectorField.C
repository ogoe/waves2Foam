/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "gabcVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "crossVersionCompatibility.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

gabcVelocityFvPatchVectorField::
gabcVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(p, iF),
    phiName_("phi"),
    isZeroVorticity_(false),
    phip_(p.size(), 0),
    timeIndex_(-1)
{
    refValue() = vector::zero;
    refGrad() = vector::zero;
    valueFraction() = symmTensor::zero;
}


gabcVelocityFvPatchVectorField::
gabcVelocityFvPatchVectorField
(
    const gabcVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    isZeroVorticity_(ptf.isZeroVorticity_),
    phip_(ptf.phip_),
    timeIndex_(ptf.timeIndex_)
{
}


gabcVelocityFvPatchVectorField::
gabcVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    isZeroVorticity_(Switch(dict.lookup("setZeroVorticity"))),
    phip_(p.size(), 0),
    timeIndex_(-1)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

    refValue() = vector::zero;
    refGrad() = vector::zero;
    valueFraction() = symmTensor::zero;
}


gabcVelocityFvPatchVectorField::
gabcVelocityFvPatchVectorField
(
    const gabcVelocityFvPatchVectorField& pivpvf
)
:
    directionMixedFvPatchVectorField(pivpvf),
    phiName_(pivpvf.phiName_),
    isZeroVorticity_(pivpvf.isZeroVorticity_),
    phip_(pivpvf.phip_),
    timeIndex_(pivpvf.timeIndex_)
{}


gabcVelocityFvPatchVectorField::
gabcVelocityFvPatchVectorField
(
    const gabcVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(pivpvf, iF),
    phiName_(pivpvf.phiName_),
    isZeroVorticity_(pivpvf.isZeroVorticity_),
    phip_(pivpvf.phip_),
    timeIndex_(pivpvf.timeIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void gabcVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    directionMixedFvPatchVectorField::autoMap(m);
}


void gabcVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    directionMixedFvPatchVectorField::rmap(ptf, addr);
}


void gabcVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = this->patch().boundaryMesh().mesh();

    // Get the face flux in the first iteration to align with the face flux
    // used in the momentum equation.
    // At first, simply phiName_ from every time step was used, but that caused
    // wiggles at the boundaries.
    // NJ 2018-11-03
    if (timeIndex_ < mesh.thisDb().time().timeIndex())
    {
    	timeIndex_ = mesh.thisDb().time().timeIndex();

    	const scalarField temp = this->patch().boundaryMesh().mesh().objectRegistry
    	    ::lookupObject<surfaceScalarField>(phiName_).boundaryField()[this->patch().index()];

    	phip_.setSize(temp.size());
    	phip_ = temp;
    }

    // Get the VOF-field on the boundary
    scalarField alphap
    (
        this->patch().boundaryMesh().mesh().objectRegistry
        ::lookupObject<volScalarField>(Foam::waves2Foam::aName())
        .boundaryField()[this->patch().index()]
    );
    alphap = Foam::min(Foam::max(alphap, 0.0), 1.0);

    scalarField scaling(Foam::pos(alphap - Foam::waves2Foam::alphaInterface()));

    // Get the normal vector
    vectorField nf = this->patch().nf();
    nf /= Foam::mag(nf);

    // Get the magnitude of the face normals (face areas)
    scalarField magSf = Foam::mag(this->patch().Sf());

    valueFraction() = sqr(nf);
    refValue() = nf*phip_/magSf;

    if (isZeroVorticity_)
    {
        // Get the vertical direction
        const uniformDimensionedVectorField& G =
            mesh.thisDb().lookupObject<uniformDimensionedVectorField>("g");

        vector vertDir(G.value()/Foam::mag(G.value()));

        // Construct the helper function and calculate the gradient
        gabcHelperFunctions gabc(mesh, vertDir, this->patch().index());
        scalarField vel(phip_/magSf*scaling);

        tmp<vectorField> gradient(gabc.normalGradient(vel, alphap));

        refGrad() = scaling*gradient();
    }
    else
    {
        refGrad() = vector::zero;
    }

    directionMixedFvPatchVectorField::updateCoeffs();
    directionMixedFvPatchVectorField::evaluate();
}


void gabcVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    if (phiName_ != "phi")
    {
        os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    }
    os.writeKeyword("setZeroVorticity") << isZeroVorticity_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void gabcVelocityFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    vectorField normalValue = transform(valueFraction(), refValue());
    vectorField transformGradValue = transform(I - valueFraction(), pvf);
    fvPatchField<vector>::operator=(normalValue + transformGradValue);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    gabcVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
