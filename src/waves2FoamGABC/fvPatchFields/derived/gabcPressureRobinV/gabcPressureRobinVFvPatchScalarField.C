/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "gabcPressureRobinVFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "crossVersionCompatibility.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

gabcPressureRobinVFvPatchScalarField::gabcPressureRobinVFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    robinVFvPatchField<scalar>(p, iF),

    waveProps_
    (
        waveTheories::waveTheory::New
        (
            this->patch().name(),
#if EXTBRANCH==1
            this->dimensionedInternalField().mesh()
#elif OFPLUSBRANCH==1
    #if OFVERSION < 1812
            this->dimensionedInternalField().mesh()
    #else
            this->internalField().mesh()
    #endif
#else
    #if OFVERSION<400
            this->dimensionedInternalField().mesh()
    #else
            this->internalField().mesh()
    #endif
#endif
        )
    ),

    shapeFunc_
    (
        Foam::celerityShapeFunctions::New
        (
#if EXTBRANCH==1
            this->dimensionedInternalField().mesh()
#elif OFPLUSBRANCH==1
    #if OFVERSION < 1812
            this->dimensionedInternalField().mesh()
    #else
            this->internalField().mesh()
    #endif
#else
    #if OFVERSION<400
            this->dimensionedInternalField().mesh()
    #else
            this->internalField().mesh()
    #endif
#endif
            ,
            this->patch().name(),
            waveProps_->patchDict(),
            waveProps_->returnDir(),
            waveProps_->seaLevel()
        )
    )
{
    checkPressure();
}


gabcPressureRobinVFvPatchScalarField::gabcPressureRobinVFvPatchScalarField
(
    const gabcPressureRobinVFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    robinVFvPatchField<scalar>(ptf, p, iF, mapper),

    waveProps_(ptf.waveProps_),

    shapeFunc_(ptf.shapeFunc_)
{
    checkPressure();
}


gabcPressureRobinVFvPatchScalarField::gabcPressureRobinVFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    robinVFvPatchField<scalar>(p, iF, dict),

    waveProps_
    (
        waveTheories::waveTheory::New
        (
           this->patch().name(),
#if EXTBRANCH==1
           this->dimensionedInternalField().mesh()
#elif OFPLUSBRANCH==1
    #if OFVERSION < 1812
            this->dimensionedInternalField().mesh()
    #else
            this->internalField().mesh()
    #endif
#else
    #if OFVERSION<400
            this->dimensionedInternalField().mesh()
    #else
            this->internalField().mesh()
    #endif
#endif
        )
    ),

    shapeFunc_
    (
        Foam::celerityShapeFunctions::New
        (
#if EXTBRANCH==1
            this->dimensionedInternalField().mesh()
#elif OFPLUSBRANCH==1
    #if OFVERSION < 1812
            this->dimensionedInternalField().mesh()
    #else
            this->internalField().mesh()
    #endif
#else
    #if OFVERSION<400
            this->dimensionedInternalField().mesh()
    #else
            this->internalField().mesh()
    #endif
#endif
            ,
            this->patch().name(),
            waveProps_->patchDict(),
            waveProps_->returnDir(),
            waveProps_->seaLevel()
        )
    )
{
    checkPressure();

    evaluate();
}


gabcPressureRobinVFvPatchScalarField::gabcPressureRobinVFvPatchScalarField
(
    const gabcPressureRobinVFvPatchScalarField& ptf
)
:
    robinVFvPatchField<scalar>(ptf),

    waveProps_(ptf.waveProps_),

    shapeFunc_(ptf.shapeFunc_)
{
    checkPressure();
}


gabcPressureRobinVFvPatchScalarField::gabcPressureRobinVFvPatchScalarField
(
    const gabcPressureRobinVFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    robinVFvPatchField<scalar>(ptf, iF),

    waveProps_(ptf.waveProps_),

    shapeFunc_(ptf.shapeFunc_)
{
    checkPressure();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void gabcPressureRobinVFvPatchScalarField::checkPressure() const
{
    if (!waveProps_->implementPressure())
    {
        FatalErrorIn("void gabcPressureRobinVFvPatchScalarField::checkPressure()")
            << "The chosen wave theory does not implement the pressure."
            << endl << endl << exit(FatalError);
    }
}


void gabcPressureRobinVFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    Info << "AUTOMAP NOT CODED!" << endl;
}


void gabcPressureRobinVFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
//    robinVFvPatchField<scalar>::rmap(ptf, addr);

    Info << "RMAP NOT CODED!" << endl;

    const gabcPressureRobinVFvPatchScalarField& mptf =
        refCast<const gabcPressureRobinVFvPatchScalarField >(ptf);
}


void gabcPressureRobinVFvPatchScalarField::updateCoeffs()
{
    // Get a const reference to the computational mesh
    const fvMesh& mesh = this->patch().boundaryMesh().mesh();
    const scalar t = mesh.thisDb().time().value();

#if OFPLUSBRANCH==1
    #if OFVERSION<1812
        scalar g(Foam::mag(uniformDimensionedVectorField
            (
                mesh.thisDb().lookupObject<uniformDimensionedVectorField>("g")
            ).value()));
    #else
        scalar g(Foam::mag((meshObjects::gravity::New(mesh.thisDb().time()).value())));
    #endif
#endif

    scalarField scaling(this->sourceValue_.size(), 1.0);

    // Only update the coefficients, if the field "Dp" is known. This is a field
    // in the pressure-correction loop.
    if (mesh.foundObject<surfaceScalarField>(Foam::waves2Foam::rAUfName()))
    {
    	// Store the index of the current patch
    	label patchID = this->patch().index();

        // Get the off-diagonal and explicit fluxes
        scalarField phiHbyA
        (
   			mesh.objectRegistry
   			::lookupObject<surfaceScalarField>("phiHbyA").boundaryField()[patchID]
        );

        scalarField phi
        (
            mesh.objectRegistry::lookupObject<surfaceScalarField>("phi")
            .boundaryField()[patchID]
        );

        // Get the inverse of the diagonal of the momentum equation
        scalarField rAUf
        (
            this->patch().boundaryMesh().mesh().objectRegistry
         	::lookupObject<surfaceScalarField>(Foam::waves2Foam::rAUfName())
         	 .boundaryField()[patchID]
        );

        // Get the density
        scalarField rho
        (
   			this->patch().boundaryMesh().mesh().objectRegistry
   			::lookupObject<volScalarField>("rho").boundaryField()[patchID]
        );

        // Get the void fraction field
        scalarField alphaGABC
        (
   			this->patch().boundaryMesh().mesh().objectRegistry
   			::lookupObject<volScalarField>(Foam::waves2Foam::aName()).boundaryField()[patchID]
        );

        // Limit the field of the void fraction field to between [0, 1]
        alphaGABC = Foam::min(Foam::max(alphaGABC, 0.0), 1.0);
        scaling = Foam::pos(alphaGABC - Foam::waves2Foam::alphaInterface());

    	// Get the magnitude of the face normals (face areas)
    	scalarField magSf = Foam::mag(this->patch().Sf());

    	// Evaluate the propagation speed
        tmp<scalarField> trhoC = shapeFunc_->rhoCelerity
            (this->patch(), rho, alphaGABC, g);
        const scalarField& rhoC = trhoC();

        // Prepare the target values
        const vectorField& cf = this->patch().Cf();
        scalarField sourceField(this->sourceValue_.size(), 0.0);

        scalarField targetP(this->sourceValue_.size(), 0.0);
        vectorField targetU(this->sourceValue_.size(), vector::zero);

        scalar rhoW = waveProps_->rhoWater();

        forAll (sourceField, facei)
        {
        	// Dummy limit, because scaling is binary (yet double)
        	if (0.5 < scaling[facei])
        	{
                targetP[facei] = waveProps_->pExcess(cf[facei], t)*rho[facei]/rhoW;
                targetU[facei] = waveProps_->U(cf[facei], t);
        	}
        }

    	vectorField nf = this->patch().nf();
        sourceField = -targetP + rhoC*(targetU & nf);

        // Evaluate the source term
    	this->sourceValue_ = pTraits<scalar>::one*
    		(rhoC/Foam::mag(magSf)*phiHbyA - sourceField);

    	// Evaluate the diffusion coefficient, which a non-unity value below
    	// the free surface
    	this->gradCoeffs_ = scaling*rhoC*rAUf + (1 - scaling)*rAUf;

        // There is only a fixed value below the free surface
        this->fixedCoeffs_ = scaling;

        // Source field is given below the free surface following wave generation
        // and as a fluxPressure condition above
        this->sourceValue_ *= scaling;
        this->sourceValue_ += (1 - scaling)*(phiHbyA - phi)/Foam::mag(magSf);
    }
    else
    {
        this->fixedCoeffs_ = 0.0;
        this->gradCoeffs_  = 1.0;
        this->sourceValue_ = 0.0;
    }
}


void gabcPressureRobinVFvPatchScalarField::write(Ostream& os) const
{
    robinVFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    gabcPressureRobinVFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
