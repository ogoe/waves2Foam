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
    Foam::jensenJacobsenChristensen2014

Description


SourceFiles
    jensenJacobsenChristensen2014.C

Author
    Niels G. Jacobsen, Deltares.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef jensenJacobsenChristensen2014_H
#define jensenJacobsenChristensen2014_H

#include "wavesPorosityModel.H"
#include "jjc2014Zones.H"
//#include "porosityZones.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class jensenJacobsenChristensen2014 Declaration
\*---------------------------------------------------------------------------*/

class jensenJacobsenChristensen2014
:
    public wavesPorosityModel
{

protected:

    // Protected data
//        porosityZones pZones_;
        jjc2014Zones pZones_;

    // Protected member functions


public:

    //- Runtime type information
    TypeName("jensenJacobsenChristensen2014");

    // Constructors

        //- from components
        jensenJacobsenChristensen2014
        (
            const fvMesh& mesh
        );


    // Destructor

        ~jensenJacobsenChristensen2014()
        {}


    // Member Functions

        //- mirror fvm::ddt with porosity
        virtual tmp<fvMatrix<vector> > ddt
        (
            GeometricField<vector, fvPatchField, volMesh>&
        );

        //- mirror fvm::ddt with porosity
        virtual tmp<fvMatrix<vector> > ddt
        (
            const geometricOneField&,
            GeometricField<vector, fvPatchField, volMesh>&
        );

        //- mirror fvm::ddt with porosity
        virtual tmp<fvMatrix<vector> > ddt
        (
            const dimensionedScalar&,
            GeometricField<vector, fvPatchField, volMesh>&
        );

        //- mirror fvm::ddt with porosity
        virtual tmp<fvMatrix<vector> > ddt
        (
            const volScalarField&,
            GeometricField<vector, fvPatchField, volMesh>&
        );

        //- Update the porosity field
        virtual void updatePorosity();

        //- Return a volScalarField with the porosities
//        virtual tmp<volScalarField> porosity() const;
        const volScalarField& porosity() const;

        //- Add the viscous and inertial resistance force contribution
        //  to the momentum equation
        virtual void addResistance(fvVectorMatrix& UEqn) const;

        //- Add the viscous and inertial resistance force contribution
        //  to the tensorial diagonal
        virtual void addResistance
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU
        ) const;



};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
