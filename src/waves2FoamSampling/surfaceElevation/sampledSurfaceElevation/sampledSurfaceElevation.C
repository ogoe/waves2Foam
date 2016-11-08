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

\*---------------------------------------------------------------------------*/

#include "sampledSurfaceElevation.H"
#include "dictionary.H"

#if EXTBRANCH==1
    #if 310<OFVERSION
        #include "foamTime.H"
    #else
        #include "Time.H"
    #endif
#elif OFPLUSBRANCH==1
    #include "Time.H"
#else
    #include "Time.H"
#endif

#include "volFields.H"
#include "ListListOps.H"
#include "SortableList.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledSurfaceElevation, 0);
}

bool Foam::sampledSurfaceElevation::verbose_ = false;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


bool Foam::sampledSurfaceElevation::checkFieldTypes()
{
    wordList fieldTypes(fieldNames_.size());

    // check files for a particular time
    if (loadFromFiles_)
    {
        forAll (fieldNames_, fieldi)
        {
            IOobject io
            (
                fieldNames_[fieldi],
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );

#if OFPLUSBRANCH == 1
    #if OFVERSION<1606
            if (io.headerOk())
            {
                fieldTypes[fieldi] = io.headerClassName();
            }
            else
            {
                fieldTypes[fieldi] = "(notFound)";
            }
    #else
            // Circumventing the check and only check for scalarField, since 
            // we know that we do not need to check for anything else.
            if (io.typeHeaderOk<volScalarField>(true))
            {
                fieldTypes[fieldi] = io.headerClassName();
            }
            else
            {
                fieldTypes[fieldi] = "(notFound)";
            }
    #endif
#else
            if (io.headerOk())
            {
                fieldTypes[fieldi] = io.headerClassName();
            }
            else
            {
                fieldTypes[fieldi] = "(notFound)";
            }
#endif
        }
    }
    else
    {
        // check objectRegistry
        forAll (fieldNames_, fieldi)
        {
            objectRegistry::const_iterator iter =
                mesh_.find(fieldNames_[fieldi]);

            if (iter != mesh_.objectRegistry::end())
            {
                fieldTypes[fieldi] = iter()->type();
            }
            else
            {
                fieldTypes[fieldi] = "(notFound)";
            }
        }
    }


    label nFields = 0;

    // classify fieldTypes
    nFields += grep(scalarFields_, fieldTypes);
//    nFields += grep(vectorFields_, fieldTypes);
//    nFields += grep(sphericalTensorFields_, fieldTypes);
//    nFields += grep(symmTensorFields_, fieldTypes);
//    nFields += grep(tensorFields_, fieldTypes);

    if (Pstream::master())
    {
        if (debug)
        {
            Pout<< "timeName = " << mesh_.time().timeName() << nl
                << "scalarFields    " << scalarFields_ << nl;
//                << "vectorFields    " << vectorFields_ << nl
//                << "sphTensorFields " << sphericalTensorFields_ << nl
//                << "symTensorFields " << symmTensorFields_ <<nl
//                << "tensorFields    " << tensorFields_ <<nl;
        }

//        if (nFields > 0)
//        {
//            if (debug)
//            {
//                Pout<< "Creating directory "
//                    << outputPath_/mesh_.time().timeName()
//                    << nl << endl;
//            }
//
//            mkDir(outputPath_/mesh_.time().timeName());
//        }
    }

    return nFields > 0;
}


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

    const PtrList<sampledSet>& sampledSets = *this;

    forAll (sampledSets, seti)
    {
        const sampledSet& samplePts = sampledSets[seti];

        // Collect data from all processors
        List<List<point> > gatheredPts(Pstream::nProcs());
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
            ListListOps::combine<List<point> >
            (
                gatheredPts, accessOp<List<point> >()
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

        // Sort curveDist and use to fill masterSamplePts
        SortableList<scalar> sortedDist(allCurveDist);
        indexSets[seti] = sortedDist.indices();


        // The constructor for coordSet has changed as of version 2.0.
        // This is taken care of using these pre-processor statements.
#if EXTBRANCH==1
        // Get reference point (note: only master has all points)
        point refPt;

        if (allPts.size())
        {
            refPt = samplePts.getRefPoint(allPts);
        }
        else
        {
            refPt = vector::zero;
        }

        masterSampledSets.set
        (
                seti,
                new coordSet
                (
                        samplePts.name(),
                        samplePts.axis(),
                        UIndirectList<point>(allPts, indexSets[seti]),
                        refPt
                )
        );
#elif OFPLUSBRANCH==1
        masterSampledSets.set
        (
            seti,
            new coordSet
            (
                samplePts.name(),
                samplePts.axis(),
                List<point>(UIndirectList<point>(allPts, indexSets[seti])),
                allCurveDist
            )
        );
#else
    #if OFVERSION < 200
        // Get reference point (note: only master has all points)
        point refPt;

        if (allPts.size())
        {
            refPt = samplePts.getRefPoint(allPts);
        }
        else
        {
            refPt = vector::zero;
        }

        masterSampledSets.set
        (
                seti,
                new coordSet
                (
                        samplePts.name(),
                        samplePts.axis(),
                        UIndirectList<point>(allPts, indexSets[seti]),
                        refPt
                )
        );
    #else
        masterSampledSets.set
        (
                seti,
                new coordSet
                (
                        samplePts.name(),
                        samplePts.axis(),
                        List<point>(UIndirectList<point>(allPts, indexSets[seti])),
                        allCurveDist
                )
        );
    #endif
#endif
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSurfaceElevation::sampledSurfaceElevation
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    PtrList<sampledSet>(),
    name_(name),
    mesh_(refCast<const fvMesh>(obr)),
    loadFromFiles_(loadFromFiles),
    outputPath_(fileName::null),
#if EXTBRANCH==1
    searchEngine_(mesh_),
#elif OFPLUSBRANCH==1
    searchEngine_(mesh_),
#else
    #if OFVERSION<210
        searchEngine_(mesh_, true),
    #else
        searchEngine_(mesh_),
    #endif
#endif

    fieldNames_(),
    interpolationScheme_(word::null),
    writeFormat_(word::null),
    surfaceElevationFilePtr_( NULL )
{
    startTime_ = dict.lookupOrDefault<scalar>("samplingStartTime", 0.0);
    nextSampleTime_ = startTime_;
    surfaceSampleDeltaT_ =
        dict.lookupOrDefault<scalar>("surfaceSampleDeltaT", -1);

    if (Pstream::parRun())
    {
        outputPath_ = mesh_.time().path()/".."/name_;
    }
    else
    {
        outputPath_ = mesh_.time().path()/name_;
    }
    if (mesh_.name() != fvMesh::defaultRegion)
    {
        outputPath_ = outputPath_/mesh_.name();
    }

//    mkDir(outputPath_ + "/" + mesh_.time().timeName() );

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


void Foam::sampledSurfaceElevation::execute()
{
    // Do nothing - only valid on write
}


void Foam::sampledSurfaceElevation::end()
{
    // Do nothing - only valid on write
}


void Foam::sampledSurfaceElevation::write()
{
    if (size() && checkFieldTypes())
    {
        sampleIntegrateAndWrite(scalarFields_);
    }
}

bool Foam::sampledSurfaceElevation::performAction()
{
    if (surfaceSampleDeltaT_ <= 10*SMALL)
    {
        // This line is needed to update the locations of the interpolation
        // lines.
    	// Note, that performAction() is only called in case the upstream
    	// time controls passes, i.e. a given timeIndex or at outputTime().
        if (mesh_.moving())
        {
    	    this->correct();
        }

        return mesh_.time().value() >= startTime_;
    }
    else
    {
        if (mesh_.time().value() < nextSampleTime_)
        {
            return false;
        }
        else
        {
            while (mesh_.time().value() > nextSampleTime_)
            {
                nextSampleTime_ += surfaceSampleDeltaT_;
            }

            // This line is needed to update the locations of the interpolation
            // lines.
            if (mesh_.moving())
            {
        	    this->correct();
            }

            return true;
        }
    }
}


void Foam::sampledSurfaceElevation::sampleIntegrateAndWrite
(
    fieldGroup<scalar>& fields
)
{
    if (fields.size() && performAction())
    {
        scalarField result(0);
        sampleAndIntegrate(scalarFields_, result);

        if (Pstream::master())
        {
            // create file if not already there, notice: this shall be
            // done on master node only
            if (surfaceElevationFilePtr_.empty())
            {
                mkDir( outputPath_ + "/" + mesh_.time().timeName() );
                surfaceElevationFilePtr_.reset
                (
                    new OFstream
                    (
                        outputPath_ + "/" + mesh_.time().timeName()
                      + "/surfaceElevation.dat"
                    )
                );

                // write header
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
                       "void Foam::sampledSurfaceElevation::sampleIntegrateAndWrite( ... )"
                    )
                    << "Output file could not be opened in " << outputPath_
                    << "/" << mesh_.time().timeName() << endl << endl
                    << exit(FatalError);
                }
            }

            if (surfaceElevationFilePtr_.valid())
            {
                surfaceElevationFilePtr_() << mesh_.time().value();

                forAll (result, seti)
                {
                    surfaceElevationFilePtr_() << tab << result[seti];
                }

                surfaceElevationFilePtr_() << endl;
            }
        }
    }
}


void Foam::sampledSurfaceElevation::sampleAndIntegrate
(
    fieldGroup<scalar>& fields,
    Field<scalar>& result
)
{
    result.setSize(0);

    if (fields.size())
    {
        bool interpolate = interpolationScheme_ != "cell";

        // Create or use existing writer
        if (fields.formatter.empty())
        {
            fields.formatter = writer<scalar>::New(writeFormat_);
        }

        // Storage for interpolated values
        PtrList<volFieldSampler<scalar> > sampledFields(fields.size());

        forAll (fields, fieldi)
        {
            if (Pstream::master() && verbose_)
            {
                Pout<< "sampledSets::sampleAndWrite: "
                    << fields[fieldi] << endl;
            }

            if (loadFromFiles_)
            {
                GeometricField<scalar, fvPatchField, volMesh> vf
                (
                    IOobject
                    (
                        fields[fieldi],
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh_
                );

                if (interpolate)
                {
                    sampledFields.set
                    (
                        fieldi,
                        new volFieldSampler<scalar>
                        (
                            interpolationScheme_,
                            vf,
                            *this
                        )
                    );
                }
                else
                {
                    sampledFields.set
                    (
                        fieldi,
                        new volFieldSampler<scalar>(vf, *this)
                    );
                }
            }
            else
            {
                if (interpolate)
                {
                    sampledFields.set
                    (
                        fieldi,
                        new volFieldSampler<scalar>
                        (
                            interpolationScheme_,
                            mesh_.lookupObject
                            <GeometricField<scalar, fvPatchField, volMesh> >
                            (fields[fieldi]),
                            *this
                        )
                    );
                }
                else
                {
                    sampledFields.set
                    (
                        fieldi,
                        new volFieldSampler<scalar>
                        (
                            mesh_.lookupObject
                            <GeometricField<scalar, fvPatchField, volMesh> >
                            (fields[fieldi]),
                            *this
                        )
                    );
                }
            }
        }

        // Combine sampled fields from processors.
        // Note: only master results are valid
        PtrList<volFieldSampler<scalar> > masterFields(sampledFields.size());
        combineSampledValues(sampledFields, indexSets_, masterFields);

        result.setSize(masterSampledSets_.size(), 0.0);

        if (Pstream::master())
        {
            forAll (masterSampledSets_, seti)
            {
                const coordSet & cs( masterSampledSets_[seti] );

                List< const Field<scalar>*> valueSets(masterFields.size());
                valueSets[0] = &masterFields[0][seti];

                List<const Field<scalar>*> columns(valueSets.size());
                columns[0] = valueSets[0];

                const Field<scalar>& alpha = *columns[0];
                scalar tolerance(0.0001);

                if (alpha.size() < 2)
                {
                    result[seti] = -GREAT;
                }
                else if
                (
                    (alpha[0] < tolerance && alpha[alpha.size()-1] < tolerance)
                    ||
                    (
                        alpha[0] > 1.0 - tolerance &&
                        alpha[alpha.size()-1] > 1.0 - tolerance
                    )
                )
                {
                    result[seti] = -GREAT;
                }
                else
                {
                    scalar value(0);
                    scalar minScalarCoord(cs.scalarCoord(0));

                    for (int pointi=0; pointi < alpha.size() - 1; pointi++)
                    {
                        value +=
                            (
                                cs.scalarCoord(pointi + 1)
                              - cs.scalarCoord(pointi)
                             )
                             *( alpha[pointi + 1] + alpha[pointi] );

                        minScalarCoord =
                            Foam::min
                            (
                                minScalarCoord,
                                cs.scalarCoord(pointi + 1)
                            );
                    }
                    value *= 0.5;

                    result[seti] = value + minScalarCoord;
                }
            }
        }
    }
}


void Foam::sampledSurfaceElevation::read(const dictionary& dict)
{
    dict_ = dict;

    fieldNames_ = wordList(dict_.lookup("fields"));

    interpolationScheme_ = "cell";
    dict_.readIfPresent("interpolationScheme", interpolationScheme_);

    dict_.lookup("setFormat") >> writeFormat_;

    scalarFields_.clear();

    PtrList<sampledSet> newList
    (
        dict_.lookup("sets"),
        sampledSet::iNew(mesh_, searchEngine_)
    );
    transfer(newList);
    combineSampledSets(masterSampledSets_, indexSets_);

    if (Pstream::master() && debug)
    {
        Pout<< "sample fields:" << fieldNames_ << nl
            << "sample sets:" << nl << "(" << nl;

        forAll (*this, si)
        {
            Pout << "  " << operator[](si) << endl;
        }
        Pout << ")" << endl;
    }
}


void Foam::sampledSurfaceElevation::correct()
{
    // Reset interpolation
	// These two lines make the moving mesh algorithms crash
	// (tested: velocityLaplacian)
	// NGJ: 16/03/2015
//    pointMesh::Delete(mesh_);
//    volPointInterpolation::Delete(mesh_);

    searchEngine_.correct();

    // A quick test has shown that this takes a lot of time on moving meshes
    // Potentially return to improve - if possible.
    // NGJ: 16/03/2015.
    PtrList<sampledSet> newList
    (
        dict_.lookup("sets"),
        sampledSet::iNew(mesh_, searchEngine_)
    );

    transfer(newList);

    combineSampledSets(masterSampledSets_, indexSets_);
}


void Foam::sampledSurfaceElevation::updateMesh(const mapPolyMesh&)
{
    correct();
}

#if EXTBRANCH==1

void Foam::sampledSurfaceElevation::movePoints(const pointField&)
{
    correct();
}
#elif OFPLUSBRANCH==1
    void Foam::sampledSurfaceElevation::movePoints(const polyMesh&)
    {
        correct();
    }

    bool Foam::sampledSurfaceElevation::timeSet()
    {
        // Do nothing
        return true;
    }
#else

    #if OFVERSION<220

    void Foam::sampledSurfaceElevation::movePoints(const pointField&)
    {
        correct();
    }

    #else

    void Foam::sampledSurfaceElevation::movePoints(const polyMesh&)
    {
        correct();
    }

        #if OFVERSION > 220
        bool Foam::sampledSurfaceElevation::timeSet()
        {
            // Do nothing
            return true;
        }
        #elif XVERSION
        bool Foam::sampledSurfaceElevation::timeSet()
        {
            // Do nothing
            return true;
        }
        #endif

    #endif
#endif


void Foam::sampledSurfaceElevation::readUpdate
(
    const polyMesh::readUpdateState state
)
{
    if (state != polyMesh::UNCHANGED)
    {
        correct();
    }
}


// ************************************************************************* //
