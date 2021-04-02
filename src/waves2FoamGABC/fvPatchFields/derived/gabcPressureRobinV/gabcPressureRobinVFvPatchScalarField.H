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

Class
    Foam::gabcPressureRobinVFvPatchScalarField

Description
    Foam::gabcPressureRobinVFvPatchScalarField

SourceFiles
    gabcPressureRobinVFvPatchScalarField.C

Author
    Niels Gjoel Jacobsen, Deltares
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef gabcPressureRobinVFvPatchScalarField_H
#define gabcPressureRobinVFvPatchScalarField_H

#include "robinVFvPatchField.H"

#include "waveTheory.H"
#include "celerityShapeFunctions.H"
//#include "gabcHelperFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class gabcPressureRobinVFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class gabcPressureRobinVFvPatchScalarField
:
    public robinVFvPatchField<scalar>
{
    // Private data

protected:

    	autoPtr<waveTheories::waveTheory> waveProps_;
    	autoPtr<Foam::celerityShapeFunctions> shapeFunc_;

        void checkPressure() const;

public:

    //- Runtime type information
    TypeName("gabcPressureRobinV");


    // Constructors

        //- Construct from patch and internal field
        gabcPressureRobinVFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        gabcPressureRobinVFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given gabcPressureRobinVFvPatchScalarField onto a new patch
        gabcPressureRobinVFvPatchScalarField
        (
            const gabcPressureRobinVFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        gabcPressureRobinVFvPatchScalarField
        (
            const gabcPressureRobinVFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar> > clone() const
        {
            return tmp<fvPatchField<scalar> >
            (
                new gabcPressureRobinVFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        gabcPressureRobinVFvPatchScalarField
        (
            const gabcPressureRobinVFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar> > clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar> >
            (
                new gabcPressureRobinVFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access
        //- Return false: this patch field is not altered by assignment
        virtual bool assignable() const
        {
            return true;
        }

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<scalar>&,
                const labelList&
            );


        // Evaluation functions

            //- Update coefficients
            virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
