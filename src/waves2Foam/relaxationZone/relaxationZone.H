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
    Foam::relaxationZone

  Description
    Base class for relaxationSchemes.
    (This is the class formerly known as "waveRelaxationOnEcstacy".)

  SourceFiles
      relaxationZone.C

  Author
  Niels Gjoel Jacobsen and Bo Terp Paulsen, The Technical University of Denmark
  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef relaxationZone_H
#define relaxationZone_H

#include "relaxationScheme.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class relaxationZone Declaration
\*---------------------------------------------------------------------------*/

class relaxationZone
{
private:

    const fvMesh& mesh_;
    volVectorField& U_;
    volScalarField& alpha_;

    wordList relaxNames_;

    List<autoPtr<relaxationSchemes::relaxationScheme> > relaxSchemePtr_;

    void resetTargetFields();

    void correctBoundaries();

protected:

    mutable volScalarField* relaxationWeightsMomentum_;

    mutable volScalarField* relaxationWeightsSource_;

    mutable volScalarField* targetAlpha_;

    mutable volVectorField* targetVelocity_;

public:
     //- Runtime type information
    TypeName("relaxationZone");

    // Constructors

        //- Construct from components
        relaxationZone
        (
            const fvMesh&,
            volVectorField&,
            volScalarField&
        );

    // Destructor

        virtual ~relaxationZone()
        {}

    // Member Functions

        virtual void correct();

        tmp<volScalarField> numericalBeach();

    // Access Functions
        const volScalarField& relaxationWeightsMomentum() const;

        const volScalarField& relaxationWeightsSource() const;

        const volScalarField& targetAlphaField() const;

        const volVectorField& targetVelocityField() const;

};


} // End namespace Foam

#endif
