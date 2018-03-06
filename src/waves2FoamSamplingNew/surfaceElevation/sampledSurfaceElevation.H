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

#ifndef sampledSurfaceElevation_H
#define sampledSurfaceElevation_H

#include "regionFunctionObject.H"
#include "sampledSet.H"
#include "volFieldsFwd.H"
#include "meshSearch.H"
#include "interpolation.H"
#include "coordSet.H"
#include "writer.H"
#include "wordReList.H"

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

class sampledSurfaceElevation
:
    public functionObjects::regionFunctionObject,
    public PtrList<sampledSet>
{
    // Private classes

        //- Class used for grouping field types
        template<class Type>
        class fieldGroup
        :
            public DynamicList<word>
        {
        public:

            //- The set formatter
            autoPtr<writer<Type>> formatter;

            //- Construct null
            fieldGroup()
            :
                DynamicList<word>(0),
                formatter(nullptr)
            {}

            //- Construct for a particular format
            fieldGroup(const word& writeFormat)
            :
                DynamicList<word>(0),
                formatter(writer<Type>::New(writeFormat))
            {}

            //- Reset format and field list
            void clear()
            {
                DynamicList<word>::clear();
                formatter.clear();
            }

            //- Assign a new formatter
            void operator=(const word& writeFormat)
            {
                formatter = writer<Type>::New(writeFormat);
            }
        };


        //- Class used for sampling volFields
        template<class Type>
        class volFieldSampler
        :
            public List<Field<Type>>
        {
            //- Name of this collection of values
            const word name_;

        public:

            //- Construct interpolating field to the sampleSets
            volFieldSampler
            (
                const word& interpolationScheme,
                const GeometricField<Type, fvPatchField, volMesh>& field,
                const PtrList<sampledSet>&
            );

            //- Construct mapping field to the sampleSets
            volFieldSampler
            (
                const GeometricField<Type, fvPatchField, volMesh>& field,
                const PtrList<sampledSet>&
            );

            //- Construct from components
            volFieldSampler
            (
                const List<Field<Type>>& values,
                const word& name
            );

            //- Return the field name
            const word& name() const
            {
                return name_;
            }
        };


    // Static data members

        //- Output verbosity
        static bool verbose_;


    // Private data

        //- Const reference to fvMesh
        const fvMesh& mesh_;

        //- Keep the dictionary to recreate sets for moving mesh cases
        dictionary dict_;

        //- Load fields from files (not from objectRegistry)
        bool loadFromFiles_;

        //- Output path
        fileName outputPath_;

        //- Mesh search engine
        meshSearch searchEngine_;


      // Read from dictonary

        //- Names of fields to sample
        wordReList fieldSelection_;

        //- Interpolation scheme to use
        word interpolationScheme_;

        //- Output format to use
        word writeFormat_;


      // Categorized scalar/vector/tensor fields
        // Retained this formatting for practical adaptation reasons (I was lazy)
        fieldGroup<scalar> scalarFields_;


      // Merging structures

        PtrList<coordSet> masterSampledSets_;
        labelListList indexSets_;

      // Output file
         autoPtr<OFstream> surfaceElevationFilePtr_;

    // Private Member Functions

        //- Combine points from all processors. Sort by curveDist and produce
        //  index list. Valid result only on master processor.
        void combineSampledSets
        (
            PtrList<coordSet>& masterSampledSets,
            labelListList& indexSets
        );

        //- Combine values from all processors.
        //  Valid result only on master processor.
        template<class T>
        void combineSampledValues
        (
            const PtrList<volFieldSampler<T>>& sampledFields,
            const labelListList& indexSets,
            PtrList<volFieldSampler<T>>& masterFields
        );

        //- Write set on master, return fileName
        template<class Type>
        fileName writeSampleFile
        (
            const coordSet& masterSampleSet,
            const PtrList<volFieldSampler<Type>>& masterFields,
            const label setI,
            const fileName& timeDir,
            const writer<Type>& formatter
        );

        template<class Type>
        void sampleAndWrite(fieldGroup<Type>& fields);


        //- Disallow default bitwise copy construct and assignment
        sampledSurfaceElevation(const sampledSurfaceElevation&) = delete;
        void operator=(const sampledSurfaceElevation&) = delete;


public:

    //- Runtime type information
    TypeName("surfaceElevation");


    // Constructors

        //- Construct from Time and dictionary
        sampledSurfaceElevation
        (
            const word& name,
            const Time& time,
            const dictionary& dict
        );

        //- Construct for given objectRegistry and dictionary
        //  allow the possibility to load fields from files
        sampledSurfaceElevation
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~sampledSurfaceElevation();


    // Member Functions

        //- Set verbosity level
        void verbose(const bool verbosity = true);

        //- Read the surfaceElevation
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Correct for mesh changes
        void correct();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&);

        //- Update for mesh point-motion
        virtual void movePoints(const polyMesh&);

        //- Update for changes of mesh due to readUpdate
        virtual void readUpdate(const polyMesh::readUpdateState state);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sampledSurfaceElevationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
