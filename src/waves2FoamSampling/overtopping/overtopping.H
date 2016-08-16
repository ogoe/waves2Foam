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
    Foam::overtopping

Description
    Computes the overtopping rates in the case of two phase flows.

SourceFiles
    overtopping.C

\*---------------------------------------------------------------------------*/

#ifndef overtopping_H
#define overtopping_H

#include "objectRegistry.H"
#include "fvMesh.H"
#include "surfaceFieldsFwd.H"
#include "volFieldsFwd.H"
#include "mapPolyMesh.H"
#include "OFstream.H"

#include "crossVersionCompatibility.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class faceSource Declaration
\*---------------------------------------------------------------------------*/

class overtopping
{



private:

    // Private member data
        word phiName_;

        word rhoPhiName_;

        scalar rho1_, rho2_, invRhoDiff_;

        const fvMesh& mesh_;

        autoPtr<OFstream> outputFilePtr_;

        autoPtr<OFstream> outputFileForcePtr_;

    // Private member functions
        bool operateOnZone( const faceZone& ) const;

        void computeAndWriteBoundary
        (
            const label&,
            const surfaceScalarField&,
            const surfaceScalarField&,
            const surfaceScalarField&,
            const surfaceVectorField&,
            vector&,
            vector&
        ) const;

        void computeAndWrite
        (
            const faceZone&,
            const surfaceScalarField&,
            const surfaceScalarField&
        );


protected:

    // Protected data


    // Protected member functions

    // Functions to be over-ridden from IOoutputFilter class

        //- Make the output file
        virtual void makeFile();

        //- Write the output file header
        virtual void writeFileHeader();

#if EXTBRANCH==1
        //- Update mesh
        virtual void updateMesh(const mapPolyMesh&);

        //- Move points
        virtual void movePoints(const pointField&);
#else
    #if OFVERSION<220
        //- Update mesh
        virtual void updateMesh(const mapPolyMesh&);

        //- Move points
        virtual void movePoints(const pointField&);
    #endif
#endif

//#if OFVERSION<220 || EXTBRANCH==1
//        //- Update mesh
//        virtual void updateMesh(const mapPolyMesh&);
//
//        //- Move points
//        virtual void movePoints(const pointField&);
//#endif

public:

    //- Run-time type information
    TypeName("overtopping");


    //- Construct from components
    overtopping
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict,
        const bool loadFromFiles = false
    );

    virtual ~overtopping();

    // Function object functions

        //- Read from dictionary
        virtual void read(const dictionary&);

        //- Calculate and write
        virtual void write();

        //- Execute
        virtual void execute() {};

        //- Execute the at the final time-loop, currently does nothing
        virtual void end() {};

#if OFPLUSBRANCH==1
        //- Update mesh
        virtual void updateMesh(const mapPolyMesh&);

        //- Move points
        virtual void movePoints(const polyMesh&);

        virtual bool timeSet();

#elif OFVERSION>=220 && EXTBRANCH==0
        //- Update mesh
        virtual void updateMesh(const mapPolyMesh&);

        //- Move points
        virtual void movePoints(const polyMesh&);
    #if OFVERSION > 220
        virtual bool timeSet();
    #elif XVERSION
        virtual bool timeSet();
    #endif
#endif

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
