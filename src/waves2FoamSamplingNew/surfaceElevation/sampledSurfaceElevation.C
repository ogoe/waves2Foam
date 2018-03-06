/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "sampledSurfaceElevation.H"
#include "dictionary.H"
#include "Time.H"
#include "volFields.H"
#include "ListListOps.H"
#include "SortableList.H"
#include "volPointInterpolation.H"
#include "mapPolyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledSurfaceElevation, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        sampledSurfaceElevation,
        dictionary
    );
}

bool Foam::sampledSurfaceElevation::verbose_ = false;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledSurfaceElevation::combineSampledSets
(
    PtrList<coordSet>& masterSampledSets,
    labelListList& indexSets
)
{
    // Combine sampleSets from processors. Sort by curveDist. Return
    // ordering in indexSets.
    // Note: only master results are valid

    masterSampledSets_.clear();
    masterSampledSets_.setSize(size());
    indexSets_.setSize(size());

    const PtrList<sampledSet>& surfaceElevation = *this;

    forAll(surfaceElevation, setI)
    {
        const sampledSet& samplePts = surfaceElevation[setI];

        // Collect data from all processors
        List<List<point>> gatheredPts(Pstream::nProcs());
        gatheredPts[Pstream::myProcNo()] = samplePts;
        Pstream::gatherList(gatheredPts);

        List<labelList> gatheredSegments(Pstream::nProcs());
        gatheredSegments[Pstream::myProcNo()] = samplePts.segments();
        Pstream::gatherList(gatheredSegments);

        List<scalarList> gatheredDist(Pstream::nProcs());
        gatheredDist[Pstream::myProcNo()] = samplePts.curveDist();
        Pstream::gatherList(gatheredDist);

        // Combine processor lists into one big list.
        List<point> allPts
        (
            ListListOps::combine<List<point>>
            (
                gatheredPts, accessOp<List<point>>()
            )
        );
        labelList allSegments
        (
            ListListOps::combine<labelList>
            (
                gatheredSegments, accessOp<labelList>()
            )
        );
        scalarList allCurveDist
        (
            ListListOps::combine<scalarList>
            (
                gatheredDist, accessOp<scalarList>()
            )
        );


        if (Pstream::master() && allCurveDist.size() == 0)
        {
            WarningInFunction
                << "Sample set " << samplePts.name()
                << " has zero points." << endl;
        }

        // Sort curveDist and use to fill masterSamplePts
        SortableList<scalar> sortedDist(allCurveDist);
        indexSets[setI] = sortedDist.indices();

        masterSampledSets.set
        (
            setI,
            new coordSet
            (
                samplePts.name(),
                samplePts.axis(),
                List<point>(UIndirectList<point>(allPts, indexSets[setI])),
                sortedDist
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSurfaceElevation::sampledSurfaceElevation
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::regionFunctionObject(name, runTime, dict),
    PtrList<sampledSet>(),
    mesh_(refCast<const fvMesh>(obr_)),
    loadFromFiles_(false),
    outputPath_(fileName::null),
    searchEngine_(mesh_),
    interpolationScheme_(word::null),
    writeFormat_(word::null),
    surfaceElevationFilePtr_(NULL)
{
    if (Pstream::parRun())
    {
        outputPath_ = mesh_.time().path()/".."/"postProcessing"/name;
    }
    else
    {
        outputPath_ = mesh_.time().path()/"postProcessing"/name;
    }
    if (mesh_.name() != fvMesh::defaultRegion)
    {
        outputPath_ = outputPath_/mesh_.name();
    }
    // Remove ".."
    outputPath_.clean();

    read(dict);
}


Foam::sampledSurfaceElevation::sampledSurfaceElevation
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjects::regionFunctionObject(name, obr, dict),
    PtrList<sampledSet>(),
    mesh_(refCast<const fvMesh>(obr)),
    loadFromFiles_(loadFromFiles),
    outputPath_(fileName::null),
    searchEngine_(mesh_),
    interpolationScheme_(word::null),
    writeFormat_(word::null),
    surfaceElevationFilePtr_(NULL)
{
    if (Pstream::parRun())
    {
        outputPath_ = mesh_.time().path()/".."/"postProcessing"/name;
    }
    else
    {
        outputPath_ = mesh_.time().path()/"postProcessing"/name;
    }
    if (mesh_.name() != fvMesh::defaultRegion)
    {
        outputPath_ = outputPath_/mesh_.name();
    }
    // Remove ".."
    outputPath_.clean();

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSurfaceElevation::~sampledSurfaceElevation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sampledSurfaceElevation::verbose(const bool verbosity)
{
    verbose_ = verbosity;
}


bool Foam::sampledSurfaceElevation::execute()
{
    return true;
}


bool Foam::sampledSurfaceElevation::write()
{
    if (size())
    {
    	// Brute force the number of fields
        const label nFields = 1;
        scalarFields_.clear();
        scalarFields_.append(Foam::waves2Foam::aName());

        if (nFields)
        {
            sampleAndWrite(scalarFields_);
        }
    }

    return true;
}


bool Foam::sampledSurfaceElevation::read(const dictionary& dict)
{
    dict_ = dict;

    bool setsFound = dict_.found("sets");

    if (setsFound)
    {
        dict.lookup("interpolationScheme") >> interpolationScheme_;
        dict.lookup("setFormat") >> writeFormat_;

        PtrList<sampledSet> newList
        (
            dict_.lookup("sets"),
            sampledSet::iNew(mesh_, searchEngine_)
        );

        transfer(newList);
        combineSampledSets(masterSampledSets_, indexSets_);

        if (this->size())
        {
            Info<< "Reading set description:" << nl;
            forAll(*this, setI)
            {
                Info<< "    " << operator[](setI).name() << nl;
            }
            Info<< endl;
        }
    }

    if (Pstream::master() && debug)
    {
        Pout<< "sample fields:" << fieldSelection_ << nl
            << "sample sets:" << nl << "(" << nl;

        forAll(*this, setI)
        {
            Pout<< "  " << operator[](setI) << endl;
        }
        Pout<< ")" << endl;
    }

    // Create the output file
    if (Pstream::master())
    {
        if (surfaceElevationFilePtr_.empty())
        {
            // Initialise the file
            mkDir(outputPath_ + "/" + mesh_.time().timeName());
            surfaceElevationFilePtr_.reset
            (
                new OFstream
                (
                    outputPath_ + "/" + mesh_.time().timeName()
                    + "/surfaceElevation.dat"
                )
            );

            // Write header
            if (surfaceElevationFilePtr_.valid())
            {
                surfaceElevationFilePtr_() << "Time";

                forAll (masterSampledSets_, seti)
                {
                    surfaceElevationFilePtr_() << tab
                            << masterSampledSets_[seti].name();
                }
                surfaceElevationFilePtr_() << endl;

                for (int coordi = 0; coordi < 3; coordi++)
                {
                    surfaceElevationFilePtr_() << -1 - coordi;

                    forAll (masterSampledSets_, seti)
                    {
                        surfaceElevationFilePtr_() << tab
                                << masterSampledSets_[seti][0].component(coordi);
                    }
                    surfaceElevationFilePtr_() << endl;
                }
            }
            else
            {
                FatalErrorIn
                (
                        "void Foam::surfaceElevation::read( ... )"
                )
                << "Output file could not be opened in " << outputPath_
                << "/" << mesh_.time().timeName() << endl << endl
                << exit(FatalError);
            }
        }
    }

    return true;
}


void Foam::sampledSurfaceElevation::correct()
{
    bool setsFound = dict_.found("sets");

    if (setsFound)
    {
        searchEngine_.correct();

        PtrList<sampledSet> newList
        (
            dict_.lookup("sets"),
            sampledSet::iNew(mesh_, searchEngine_)
        );
        transfer(newList);
        combineSampledSets(masterSampledSets_, indexSets_);
    }
}


void Foam::sampledSurfaceElevation::updateMesh(const mapPolyMesh& mpm)
{
    if (&mpm.mesh() == &mesh_)
    {
        correct();
    }
}


void Foam::sampledSurfaceElevation::movePoints(const polyMesh& mesh)
{
    if (&mesh == &mesh_)
    {
        correct();
    }
}


void Foam::sampledSurfaceElevation::readUpdate(const polyMesh::readUpdateState state)
{
    if (state != polyMesh::UNCHANGED)
    {
        correct();
    }
}


// ************************************************************************* //
