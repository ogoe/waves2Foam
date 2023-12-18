/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::surfaceElevation

Description
    Samples the surface elevation based on the VOF-field.
    Adapted from sampledSets-class in OF+1712

    Set of sets to sample.
    Call surfaceElevation.write() to sample&write files.

SourceFiles
    surfaceElevation.C

\*---------------------------------------------------------------------------*/

#ifndef overtopping_H
#define overtopping_H

#include "regionFunctionObject.H"
//#include "sampledSet.H"
#include "volFieldsFwd.H"
#include "fvMesh.H"
//#include "meshSearch.H"
//#include "interpolation.H"
//#include "coordSet.H"
//#include "writer.H"
//#include "wordReList.H"

#include "crossVersionCompatibility.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Time;
class objectRegistry;
class dictionary;
class fvMesh;

/*---------------------------------------------------------------------------*\
                         Class surfaceElevation Declaration
\*---------------------------------------------------------------------------*/

class overtopping
:
    public functionObjects::regionFunctionObject
{
    // Private data

        //- Reference to Time
        const Time& runTime_;

        //- Output path
        fileName outputPath_;

        //- Fields
        word phiName_;
        word rhoPhiName_;
        scalar rho1_, rho2_, invRhoDiff_;

        // Output file
        autoPtr<OFstream> overtoppingFilePtr_;

    // Private Member Functions
        bool operateOnZone(const faceZone&) const;

        vector calcOvertopping
        (
            const fvMesh&,
            const surfaceScalarField&,
            const surfaceScalarField&,
            const faceZone&
        ) const;

        //- Disallow default bitwise copy construct and assignment
        overtopping(const overtopping&) = delete;
        void operator=(const overtopping&) = delete;


public:

    //- Runtime type information
    TypeName("overtopping");


    // Constructors

        //- Construct from Time and dictionary
        overtopping
        (
            const word&,
            const Time&,
            const dictionary&
        );

    //- Destructor
    virtual ~overtopping();


    // Member Functions

    //- Execute, currently does nothing
    virtual bool execute();

    //- Sample and write
    virtual bool write();

    // Read and prepare header
    bool read(const dictionary& );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
