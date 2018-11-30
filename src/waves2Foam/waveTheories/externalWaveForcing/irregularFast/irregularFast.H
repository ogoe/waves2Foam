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
    Foam::waveTheories::irregularFast

Description


SourceFiles
    irregularFast.C

Author
    Niels Gjoel Jacobsen, Deltares.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef irregularFast_H
#define irregularFast_H

#include "externalWaveForcing.H"
#include "relaxationShape.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

/*---------------------------------------------------------------------------*\
                           Class irregularFast Declaration
\*---------------------------------------------------------------------------*/

class irregularFast
:
    public externalWaveForcing
{

protected:

    // Protected data
    	// Two necessary dictionaries
        const IOdictionary& waveProps_;
        const dictionary& coeffDict_;

        // Name of the relaxation zones, which contain the external forcing
        bool ignoreMeshMotion_;

        wordList relaxationNames_;
        labelList externalSourceZones_;

        // Irregular wave properties
        label N_;
        scalar h_;

        scalarField amp_;
        scalarField omega_;
        scalarField phi_;
        vectorField k_;
        scalarField K_;
        vectorField compDir_;
        scalarField period_;
        scalarField velAmp_;
        scalar Tsoft_;

        scalar seaLevel_;

        scalar rampFactor_;

        // Parameters used in the indexing of the irregular wave field
        vector eX_, eY_, eZ_;

        scalarField verticalCoordinate_;
        List<scalarField> horizontalCoordinates_;

        scalar tolerance_;

        // The primary fields to be pre-computed
        List<scalarField> Sinh_;
        List<scalarField> Cosh_;

        List<List<scalarField> > SinK_;
        List<List<scalarField> > CosK_;

        scalarField SinOmega_;
        scalarField CosOmega_;

        List<scalarField> Eta_;

        List<List<vectorField> > Ucomponent_;
        List<List<scalarField> > Wcomponent_;


    // Protected member functions
        scalar factor(const scalar&) const;

        void calcUnitVectors();

        void calcUniqueCoordinates();

        void uniqueSorting(scalarField&);

        void calcAllFields();

        void calcHyperbolicFunctions();

        void calcWaveNumberTrigonometricFunctions();

        void calcWaveFrequencyTrigonometricFunctions();

        void calcSurfaceElevation();

        void calcVelocityComponents();

        label bisectionFindIndex
        (
            const scalar&,
            const scalarField&
        ) const;

        void findIndexing(const point&, label&, label&, label&, bool) const;

        scalar interpolateSurfaceElevation
        (
        	const point&,
        	const scalar&,
        	const label&,
        	const label&
        ) const;

        scalar surfaceElevation(const label&, const label&) const;

public:

    //- Runtime type information
    TypeName("irregularFast");

    // Constructors

        //- from components
        irregularFast
        (
            IOobject,
            Time&,
            const fvMesh&
        );


    // Destructor

        ~irregularFast()
        {}


    // Member Functions
        virtual void step();


        virtual inline void close()
        {
        	// Nothing to be done
        };


        virtual scalar eta
        (
            const point&,
            const scalar&
        ) const;


//        virtual scalar ddxPd
//        (
//            const point&,
//            const scalar&,
//            const vector&
//        ) const;


        virtual vector U
        (
            const point&,
            const scalar&
        ) const;


        virtual scalar pExcess
        (
            const point&,
            const scalar&
        ) const;

        bool implementPressure() const
        {
            return false;
        };

        virtual bool writeData(Foam::Ostream&) const
        {
            return true;
        };

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
