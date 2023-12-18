/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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
    Foam::sampledSurfaceElevation

Description
    Set of sets to sample.

    The write() method is used to sample and write files.

    Example of function object specification:

    \verbatim
    surfaces
    {
        type    sets;
        libs    (sampling);

        // Write at same frequency as fields
        writeControl    writeTime;
        writeInterval   1;

        // Fields to be sampled
        fields          (p U);

        // Scheme to obtain values
        interpolationScheme cellPoint;

        // Output format
        setFormat   vtk;

        formatOptions
        {
            vtk
            {
                precision  10;
            }
        }

        sets
        {
            Uref
            {
                type    cloud;
                axis    xyz;
                points  ((-0.0508 0.0508 0.01));
            }
        }
    }
    \endverbatim

    Entries:
    \table
        Property | Description                              | Required | Default
        type     | Type-name: sets                          | yes |
        sets     | Dictionary or list of sample sets        | expected |
        fields   | word/regex list of fields to sample      | yes |
        interpolationScheme | scheme to obtain values       | no  | cellPoint
        setFormat | output format                           | yes |
        sampleOnExecute | Sample (store) on execution as well | no | false
        formatOptions | dictionary of format options        | no  |
    \endtable

    Additional per-set entries:
    \table
        Property | Description                              | Required | Default
        store    | Store surface/fields on registry         | no  |
        setFormat | output format                           | no  |
        formatOptions | dictionary of format options        | no  |
    \endtable

Note
    Special setFormat \c probes can be used to output ensemble results
    in a format similar to the probes function object.

SourceFiles
    sampledSurfaceElevation.C
    sampledSurfaceElevationImpl.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_sampledSurfaceElevation_H
#define Foam_sampledSurfaceElevation_H

#include "fvMeshFunctionObject.H"
#include "HashPtrTable.H"
#include "OFstream.H"
#include "sampledSet.H"
#include "meshSearch.H"
#include "coordSet.H"
#include "coordSetWriter.H"
#include "volFieldsFwd.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Time;
class dictionary;
class globalIndex;

/*---------------------------------------------------------------------------*\
                         Class sampledSurfaceElevation Declaration
\*---------------------------------------------------------------------------*/

class sampledSurfaceElevation
:
    public functionObjects::fvMeshFunctionObject,
    public PtrList<sampledSet>
{
    // Static Data Members

        //- Local control for sampling actions
        enum sampleActionType : unsigned
        {
            ACTION_NONE  = 0,
            ACTION_WRITE = 0x1,
            ACTION_STORE = 0x2,
            ACTION_ALL = 0xF
        };

    // Private Data

        //- Keep the dictionary to recreate sets for moving mesh cases
        dictionary dict_;

        //- Load fields from files (not from objectRegistry)
        bool loadFromFiles_;

        //- Output verbosity
        bool verbose_;

        //- Perform sample/store actions on execute as well
        bool onExecute_;

        //- Correct meshSearch and update sets
        bool needsCorrect_;

        //- Output in raw format like 'probes' does
        bool writeAsProbes_;

        //- Output path
        fileName outputPath_;

        //- Mesh search engine
        meshSearch searchEngine_;


    // Read from dictionary

        //- Names of fields to sample
        wordRes fieldSelection_;

        //- Interpolation/sample scheme to obtain values at the points
        word samplePointScheme_;

        //- Output format to use
        word writeFormat_;

        //- Dictionary containing writer options
        dictionary writeFormatOptions_;

        // Output file
        autoPtr<OFstream> surfaceElevationFilePtr_;

    // Output control

        //- Current list of field names selected for sampling
        DynamicList<word> selectedFieldNames_;

        //- The coordSet writers (one per sampled set)
        PtrList<coordSetWriter> writers_;

        //- Current open files (non-empty on master only) when writing as probes
        HashPtrTable<OFstream> probeFilePtrs_;


    // Merging

        //- Gathered coordSet. (Content only meaningful on master)
        PtrList<coordSet> gatheredSets_;

        //- Indexed sort order for gathered coordSet.
        //  (Content only meaningful on master)
        List<labelList> gatheredSorting_;

        //- The globalIndex for gathering. (Content only meaningful on master)
        List<globalIndex> globalIndices_;


    // Private Member Functions

        //- Get or create new stream as required (on master)
        OFstream* createProbeFile(const word& fieldName);

        //- A new coordSet writer, with per-set formatOptions
        static autoPtr<coordSetWriter> newWriter
        (
            word writeType,
            const dictionary& formatOptions,
            const dictionary& surfDict
        );

        //- Perform sampling action with store/write
        bool performAction(unsigned request);

        //- Count selected/sampled fields
        //  Returns empty IOobjectList if loadFromFiles_ is not active
        //
        //  Adjusts selectedFieldNames_
        IOobjectList preCheckFields(unsigned request);

        //- Setup the sets (optional writers)
        void initDict(const dictionary& dict, const bool initial);

        //- Combine points from all processors.
        //- Sort by curve distance and produce index list.
        //- Valid result only on master processor.
        void gatherAllSets();

        //- Get from registry or load from disk
        template<class GeoField>
        tmp<GeoField> getOrLoadField(const word& fieldName) const;

        //- Sample and store/write a specific volume field
        template<class Type>
        void performAction(const VolumeField<Type>& field, unsigned request);

        //- Sample and write all applicable sampled fields
        //  Only uses IOobjectList when loadFromFiles_ is active
        template<class GeoField>
        void performAction(const IOobjectList& objects, unsigned request);

        //- No copy construct
        sampledSurfaceElevation(const sampledSurfaceElevation&) = delete;

        //- No copy assignment
        void operator=(const sampledSurfaceElevation&) = delete;


public:

    //- Runtime type information
    TypeName("surfaceElevation");


    // Constructors

        //- Construct from Time and dictionary
        sampledSurfaceElevation
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct for given objectRegistry and dictionary
        //  allow the possibility to load fields from files
        sampledSurfaceElevation
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~sampledSurfaceElevation() = default;


    // Member Functions

        //- Enable/disable verbose output
        //  \return old value
        bool verbose(const bool on) noexcept;

        //- Read the sampledSurfaceElevation
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

#endif

// ************************************************************************* //
