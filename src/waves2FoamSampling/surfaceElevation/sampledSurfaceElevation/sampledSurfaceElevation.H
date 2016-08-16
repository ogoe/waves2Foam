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
    Foam::sampledSurfaceElevation

Description
    Integrate along sampling lines (ought to be vertical), and output the
    surface elevation into the root case.

    Based on sampledSets.C. Modified to be compatible to the following
    versions:

    OF-1.6
    OF-1.7
    OF-2.0
    OF-2.1
    OF-2.2

SourceFiles
    sampledSurfaceElevation.C

Authors
    Niels G. Jacobsen and Sopheak Seng
    Technical University of Denmark

\*---------------------------------------------------------------------------*/

#ifndef sampledSurfaceElevation_H
#define sampledSurfaceElevation_H

#include "sampledSet.H"
#include "volFieldsFwd.H"
#include "meshSearch.H"
#include "interpolation.H"
#include "coordSet.H"
#include "writer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class objectRegistry;
class dictionary;
class fvMesh;

/*---------------------------------------------------------------------------*\
                      Class sampledSets Declaration
\*---------------------------------------------------------------------------*/

class sampledSurfaceElevation
:
    public PtrList<sampledSet>
{
    // Private classes

        //- Class used for grouping field types
        template<class Type>
        class fieldGroup
        :
            public wordList
        {
        public:

            //- Set formatter
            autoPtr<writer<Type> > formatter;

            //- Construct null
            fieldGroup()
            :
                wordList(0),
                formatter(NULL)
            {}

            void clear()
            {
                wordList::clear();
                formatter.clear();
            }
        };


        //- Class used for sampling volFields
        template <class Type>
        class volFieldSampler
        :
            public List<Field<Type> >
        {
            //- Name of this collection of values
            const word name_;

        public:

            //- Construct interpolating field to the sampledSurfaceElevation
            volFieldSampler
            (
                const word& interpolationScheme,
                const GeometricField<Type, fvPatchField, volMesh>& field,
                const PtrList<sampledSet>&
            );

            //- Construct mapping field to the sampledSurfaceElevation
            volFieldSampler
            (
                const GeometricField<Type, fvPatchField, volMesh>& field,
                const PtrList<sampledSet>&
            );

            //- Construct from components
            volFieldSampler
            (
                const List<Field<Type> >& values,
                const word& name
            );

            //- Return the field name
            const word& name() const
            {
                return name_;
            }
        };


    // Static data members

        //- output verbosity
        static bool verbose_;


    // Private data

        //- Name of this set of sets,
        //  Also used as the name of the sampledSets directory.
        word name_;

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

        scalar startTime_;
        scalar nextSampleTime_;
        scalar surfaceSampleDeltaT_;


        // Read from dictonary

            //- Names of fields to sample
            wordList fieldNames_;

            //- Interpolation scheme to use
            word interpolationScheme_;

            //- Output format to use
            word writeFormat_;


        // Categorized scalar/vector/tensor fields

            fieldGroup<scalar> scalarFields_;


        // Merging structures

            PtrList<coordSet> masterSampledSets_;
            labelListList indexSets_;

        // Output file
            autoPtr<OFstream> surfaceElevationFilePtr_;

    // Private Member Functions

        //- Classify field types, return true if nFields > 0
        bool checkFieldTypes();

        //- Find the fields in the list of the given type, return count
        template<class Type>
        label grep
        (
            fieldGroup<Type>& fieldList,
            const wordList& fieldTypes
        ) const;

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
            const PtrList<volFieldSampler<T> >& sampledFields,
            const labelListList& indexSets,
            PtrList<volFieldSampler<T> >& masterFields
        );

        void sampleIntegrateAndWrite( fieldGroup<scalar>& fields );

        void sampleAndIntegrate( fieldGroup<scalar>& fields, Field<scalar>& );

        //- Disallow default bitwise copy construct and assignment
        sampledSurfaceElevation(const sampledSurfaceElevation&);
        void operator=(const sampledSurfaceElevation&);

        bool performAction();

public:

    //- Runtime type information
    TypeName("surfaceElevation");


    // Constructors

        //- Construct for given objectRegistry and dictionary
        //  allow the possibility to load fields from files
        sampledSurfaceElevation
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    // Destructor

        virtual ~sampledSurfaceElevation();


    // Member Functions

        //- Return name of the set of probes
        virtual const word& name() const
        {
            return name_;
        }

        //- set verbosity level
        void verbose(const bool verbosity = true);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Sample and write
        virtual void write();

        //- Read the sampledSets
        virtual void read(const dictionary&);

        //- Correct for mesh changes
        void correct();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&);

        //- Update for mesh point-motion

#if EXTBRANCH==1
        virtual void movePoints(const pointField&);
#elif OFPLUSBRANCH==1
        virtual void movePoints(const polyMesh&);

        virtual bool timeSet();
#else
    #if OFVERSION<220
        virtual void movePoints(const pointField&);
    #else
        virtual void movePoints(const polyMesh&);

        #if OFVERSION>220
            virtual bool timeSet();
        #elif XVERSION
            virtual bool timeSet();
        #endif
    #endif
#endif

        //- Update for changes of mesh due to readUpdate
        virtual void readUpdate(const polyMesh::readUpdateState state);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "sampledSurfaceElevationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
