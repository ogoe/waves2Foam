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
    Foam::wavesPorosityModel

Description


SourceFiles
    wavesPorosityModel.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef wavesPorosityModel_H
#define wavesPorosityModel_H

#include "fvMesh.H"
#include "volFields.H"
#include "IOdictionary.H"

#include "autoPtr.H"
#include "runTimeSelectionTables.H"

#include "uniformDimensionedFields.H"

#include "fvCFD.H"
#include "volFieldsFwd.H"
#include "fvMatrix.H"
#include "geometricOneField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class wavesPorosityModel Declaration
\*---------------------------------------------------------------------------*/

class wavesPorosityModel
{

protected:

    // Protected data
        mutable volScalarField porosity_;

    // Protected member functions


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        wavesPorosityModel(const wavesPorosityModel&);

        //- Disallow default bitwise assignment
        void operator=(const wavesPorosityModel&);

public:

    //- Runtime type information
    TypeName("wavesPorosityModel");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            wavesPorosityModel,
            wavesPorosityModel,
            (
                const fvMesh& mesh
            ),
            (mesh)
        );

    // Constructors

        //- Construct from components
        wavesPorosityModel
        (
            const fvMesh&
        );

    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<wavesPorosityModel> New
        (
            const fvMesh&
        );

    // Destructor

        virtual ~wavesPorosityModel();

    // Member Functions
        //- The templated part had to be skipped, because it should also be a
        //  virtual function. I could not get this part working.
        //  NGJ, 29.10.2014

        //- mirror fvm::ddt with porosity
        virtual tmp<fvMatrix<vector> > ddt
        (
            GeometricField<vector, fvPatchField, volMesh>&
        ) = 0;

        //- mirror fvm::ddt with porosity
        virtual tmp<fvMatrix<vector> > ddt
        (
            const geometricOneField&,
            GeometricField<vector, fvPatchField, volMesh>&
        ) = 0;

        //- mirror fvm::ddt with porosity
        virtual tmp<fvMatrix<vector> > ddt
        (
            const dimensionedScalar&,
            GeometricField<vector, fvPatchField, volMesh>&
        ) = 0;

        //- mirror fvm::ddt with porosity
        virtual tmp<fvMatrix<vector> > ddt
        (
            const volScalarField&,
            GeometricField<vector, fvPatchField, volMesh>&
        ) = 0;

        //- Update the porosity field
        virtual void updatePorosity() = 0;

        //- Return a volScalarField with the porosities
//        virtual tmp<volScalarField> porosity() const = 0;
        virtual const volScalarField& porosity() const = 0;

        //- Add the viscous and inertial resistance force contribution
        //  to the momentum equation
        virtual void addResistance(fvVectorMatrix& UEqn) const = 0;

        //- Add the viscous and inertial resistance force contribution
        //  to the tensorial diagonal
        virtual void addResistance
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU
        ) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
