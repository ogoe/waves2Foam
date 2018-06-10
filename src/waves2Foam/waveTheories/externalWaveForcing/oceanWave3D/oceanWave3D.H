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
    Foam::waveTheories::oceanWave3D

Description


SourceFiles
    oceanWave3D.C

Author
    Niels Gjoel Jacobsen, Deltares.  All rights reserved.
    Bo Terp Paulsen, Deltares. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef oceanWave3D_H
#define oceanWave3D_H

#include "externalWaveForcing.H"
#include "relaxationShape.H"

#include "crossVersionCompatibility.H"
#include "setWaveField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

/*---------------------------------------------------------------------------*\
                           Class oceanWave3D Declaration
\*---------------------------------------------------------------------------*/

class oceanWave3D
:
    public externalWaveForcing
{

protected:

    // Protected data
    	// Two necessary dictionaries
        const IOdictionary& waveProps_;
        const dictionary& coeffDict_;

        scalar seaLevel_;

        label N_;
        scalarField startTimes_;
        scalarField endTimes_;

        scalar startTime_;
        scalar endTime_;

        Switch ramp_;
        scalar Tsoft_;

        scalar maxDT_;

        vector translateOFMesh_;

        tensor OFtoOCW_;
        tensor OCWtoOF_;


    // Protected member functions
        scalar factor(const scalar&) const;

        //- Function used, when OCW3D alone takes a time step.
        dimensionedScalar OCWTimeStep() const;

        //- Construct mapping tensors, if vertical coordinate does not match
        void mappingTensors();

        //- Update the OpenFoam intervals
        void updateIntervals();

        //- Align the time instances between OceanWave3D and OpenFoam
        void alignTimes();

        //- This function time steps OCw3D until it reaches the beginning of
        //  an OpenFoam time interval
        void timeStepOceanWave3D();

        //- Map the OceanWave3D solution into the OpenFoam domain
        void mapSolution();

        //- Update the face flux field by solving correctPhi.H
        void updatePhi();

        //- Update the time steps according to the OF-field: phi
        void updateTimeAndTimeStep();

        void takeTimeStep(bool);

        void takeTimeStep(const scalar, bool);

        void calculateKinematics();

        void writeExternal() const;


public:

    //- Runtime type information
    TypeName("oceanWave3D");

    // Constructors

        //- from components
        oceanWave3D
        (
            IOobject,
            Time&,
            const fvMesh&
        );


    // Destructor

        ~oceanWave3D()
        {}


    // Member Functions
        virtual void step();


        virtual inline void close();


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
        	Info << "WRITING OCW3D-DATA" << endl;
            return true;
        };

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
