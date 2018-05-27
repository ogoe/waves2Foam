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
    Foam::waveTheories::externalWaveForcing

Description


SourceFiles
    externalWaveForcing.C

Author
    Niels Gjoel Jacobsen, Deltares.  All rights reserved.
    Bo Terp Paulsen, Deltares, All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef externalWaveForcing_H
#define externalWaveForcing_H

#include "fvMesh.H"
#include "volFields.H"
#include "IOdictionary.H"

#include "autoPtr.H"
#include "runTimeSelectionTables.H"

#include "uniformDimensionedFields.H"
#include "crossVersionCompatibility.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

/*---------------------------------------------------------------------------*\
                         Class externalWaveForcing Declaration
\*---------------------------------------------------------------------------*/

class externalWaveForcing
:
    public regIOobject
{

protected:

    // Protected data
        Time& rT_;

        const fvMesh& mesh_;

    // Protected member functions


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        externalWaveForcing(const externalWaveForcing&);

        //- Disallow default bitwise assignment
        void operator=(const externalWaveForcing&);

public:

    //- Runtime type information
    TypeName("externalWaveForcing");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            externalWaveForcing,
            externalWaveForcing,
            (
                IOobject io,
                Time& rT,
                const fvMesh& mesh
            ),
            (io, rT, mesh)
        );

    // Constructors

        //- Construct from components
        externalWaveForcing
        (
            IOobject io,
            Time& rT,
            const fvMesh& mesh
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<externalWaveForcing> New
        (
            IOobject io,
            Time& rT,
            const fvMesh& mesh
        );

    // Destructor

        virtual ~externalWaveForcing();

    // Member Functions
        virtual void step() = 0;


        virtual void close() = 0;


        virtual scalar eta
        (
            const point&,
            const scalar&
        ) const = 0;


//        virtual scalar ddxPd
//        (
//            const point&,
//            const scalar&,
//            const vector&
//        ) const = 0;

        virtual scalar pExcess
        (
            const point&,
            const scalar&
        ) const = 0;

        virtual bool implementPressure() const = 0;

        virtual vector U
        (
            const point&,
            const scalar&
        ) const = 0;

        virtual bool writeData(Foam::Ostream&) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
