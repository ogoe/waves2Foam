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

#include "postProcessingWaves.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


bool pairSortA( std::pair<scalar, label> a, std::pair<scalar, label> b)
{
    return a.first < b.first;
}


bool pairSortB( std::pair<scalar, word> a, std::pair<scalar, word> b)
{
    return a.first < b.first;
}


namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(postProcessingWaves, 0);
defineRunTimeSelectionTable(postProcessingWaves, postProcessingWaves);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //


void postProcessingWaves::getTimeDirs
(
    const word& inputDir,
    wordList& timeDirs
)
{
    if (actionProperties_.found("timeDirs"))
    {
        wordList temp( actionProperties_.lookup("timeDirs") );
        timeDirs = temp;
    }
    else
    {
        // Get the time directories
        fileNameList fnl =
            Foam::readDir( inputDir, Foam::fileName::DIRECTORY );

        timeDirs.setSize( fnl.size() );

        forAll (timeDirs, timei)
        {
            timeDirs[timei] = fnl[timei].name();
        }
    }

    // Make a numerical sorting of the words
    // This sorting is needed to easily exclude duplicate time instances in
    // the outputted data, when read by raw-data methods
    List<std::pair<scalar, word> > timeWord( timeDirs.size() );

    forAll (timeWord, timei)
    {
        timeWord[timei].first = std::atof( timeDirs[timei].c_str() );
        timeWord[timei].second = timeDirs[timei];
    }

    std::sort( timeWord.begin(), timeWord.end(), pairSortB );

    forAll (timeDirs, timei)
    {
        timeDirs[timei] = timeWord[timei].second;
    }
}


void postProcessingWaves::writeNameDict
(
    const scalar& dt,
    const wordList& names
)
{
    IOdictionary xyz
    (
        IOobject
        (
            callName_ + "_dict",
            rT_.constant(),
            addDir_,
            rT_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    );

    // Adding deltaT information to the dictionary
    xyz.add("deltaT", dt, true);

    // Adding indexing to the dictionary
    labelList index( names.size(), 0 );

    forAll (index, indexi)
    {
        index[indexi] = indexi;
    }

    xyz.add("index", index, true);

    // Adding the names
    xyz.add("names", names, true);

    // Write the dictionary
    xyz.regIOobject::write();
}


void postProcessingWaves::writeXYZDict
(
    const scalar& dt,
    const scalarField& x,
    const scalarField& y,
    const scalarField& z
)
{
    // Open a dictionary used for outputting data
    IOdictionary xyz
    (
        IOobject
        (
            callName_ + "_dict",
            rT_.constant(),
            addDir_,
            rT_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    );

    // Adding deltaT information to the dictionary
    xyz.add("deltaT", dt, true);

    // Adding indexing to the dictionary
    labelList index( x.size(), 0 );

    forAll (index, indexi)
    {
        index[indexi] = indexi;
    }

    xyz.add("index", index, true);

    // Adding the point locations to the dictionary
    xyz.add("x", x, true);

    xyz.add("y", y, true);

    xyz.add("z", z, true);

    // Write the dictionary
    xyz.regIOobject::write();
}


void postProcessingWaves::writeIOScalarField
(
    const scalarField& field,
    const word& name
)
{
    IOField<scalar> output
    (
        IOobject
        (
            name,
            rT_.constant(),
            addDir_,
            rT_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        field
    );
    output.write();
}


void postProcessingWaves::writeIOVectorField
(
    const vectorField& field,
    const word& name
)
{
    IOField<vector> output
    (
        IOobject
        (
            name,
            rT_.constant(),
            addDir_,
            rT_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        field
    );
    output.write();
}


scalarField postProcessingWaves::readIOScalarField
(
    const word& name
)
{
    IOField<scalar> field
    (
        IOobject
        (
            name,
            rT_.constant(),
            addDir_,
            rT_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    return field;
}


vectorField postProcessingWaves::readIOVectorField
(
    const word& name
)
{
    IOField<vector> field
    (
        IOobject
        (
            name,
            rT_.constant(),
            addDir_,
            rT_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    return field;
}


word postProcessingWaves::dataType()
{
    std::stringstream ss;
    ss << callName_ << "_0";

    IOobject fileHeader
    (
        ss.str(),
        rT_.constant(),
        addDir_,
        rT_,
        IOobject::MUST_READ
    );

    // The inquiry to fileHeader.headerOk() is needed to update
    // headerClassName() from IOobject. This is weird but not looked into.
#if OFPLUSBRANCH == 1
    #if OFVERSION<1606
        if (!fileHeader.headerOk())
        {
        }
    #else
        // Dummy question to read in the class name
        fileHeader.typeHeaderOk<volScalarField>(false);
    #endif
#else
    if (!fileHeader.headerOk())
    {
    }
#endif

    if (fileHeader.headerClassName() == "scalarField")
    {
        return "scalar";
    }
    else if (fileHeader.headerClassName() == "vectorField")
    {
        return "vector";
    }
    else
    {
        return "unknown";
    }
}


scalar postProcessingWaves::readDeltaT
(
    const dictionary& timeDict
)
{
    scalar dt( readScalar(timeDict.lookup("deltaT")) );

    // Check for validity of the time step
    if (dt <= 0.0)
    {
        FatalErrorIn("scalar postProcessingWaves::readDeltaT")
            << "The time step (deltaT) given in\n    "
            <<  timeDict.name() << endl
            << "is less than or equal to zero." << endl << exit(FatalError);
    }

    return dt;
}


void postProcessingWaves::readIndices
(
    const dictionary& indexDict,
    labelList& indices
)
{
    // Getting the labelList of data set indices
    if (actionProperties_.lookupOrDefault<Switch>("allDataSets", false))
    {
        labelList tempList( indexDict.lookup("index") );
        indices.setSize( tempList.size() );
        indices = tempList;
    }
    else
    {
        labelList tempList( actionProperties_.lookup("indices") );
        indices.setSize( tempList.size() );
        indices = tempList;
    }
}


List<scalarField> postProcessingWaves::readScalarFields
(
    const labelList& indices
)
{
    List<scalarField> res( indices.size() );

    forAll (indices, indexi)
    {
        std::stringstream ss;
        ss << callName_ << "_" << indices[indexi];

        scalarField& r( res[indexi] );

        r = readIOScalarField( ss.str() );
    }

    return res;
}


List<vectorField> postProcessingWaves::readVectorFields
(
    const labelList& indices
)
{
    List<vectorField> res( indices.size() );

    forAll (indices, indexi)
    {
        std::stringstream ss;
        ss << callName_ << "_" << indices[indexi];

        vectorField& r( res[indexi] );

        r = readIOVectorField( ss.str() );
    }

    return res;
}


scalarField postProcessingWaves::equidistantTime
(
    const List<std::pair<scalar, label > >& timeLabel,
    const dictionary& dict
)
{
    scalar tmin( dict.lookupOrDefault<scalar>("tMin", timeLabel[0].first ) );
    scalar tmax
    (
        dict.lookupOrDefault<scalar>
        (
            "tMax",
            timeLabel[timeLabel.size() - 1].first
        )
    );
    scalar dt( readScalar( dict.lookup("deltaT" ) )  );

    Info << "        - Interpolation range: [" << tmin << "; " << tmax << "]"
         << endl;

    label N( static_cast<label>( ( tmax - tmin )/dt ) );

    scalarField t(N, tmin);

    forAll (t, timei)
    {
        t[timei] += timei*dt;
    }

    return t;
}


void postProcessingWaves::interpolationWeights
(
    const List<std::pair<scalar, label> >& timeLabel,
    const scalarField& t,
    scalarField& weights,
    labelList& first,
    labelList& second
)
{
    weights.setSize( t.size() );
    first.setSize( t.size() );
    second.setSize( t.size() );

    if
    (
        t[0] < timeLabel[0].first ||
        timeLabel[timeLabel.size() - 1].first < t[t.size() - 1]
    )
    {
        FatalErrorIn
        (
            "void Foam::dataProcessingTools::interpolationWeights( ... )"
        )
        << "The target times for interpolations are outside the bounds" << endl
        << "of the sampled data set." << endl << exit(FatalError);
    }

    label count(0);

    forAll (t, timei)
    {
        while
        (
            t[timei] >= timeLabel[count].first &&
            t[timei] >= timeLabel[count + 1].first
        )
        {
            count++;
        }

        scalar w
        (
            (timeLabel[count + 1].first - t[timei])
           /(timeLabel[count + 1].first - timeLabel[count].first)
        );

        weights[timei] = w;
        first[timei]   = timeLabel[count].second;
        second[timei]  = timeLabel[count + 1].second;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


postProcessingWaves::postProcessingWaves
(
    const Time& rT,
    const dictionary& actionProp,
    const word& action
)
:
    rT_(rT),

    actionProperties_( actionProp ),

    callName_( actionProperties_.lookup("callName") ),

    actionType_( action ),

    addDir_( "../postProcessedWaves2Foam/" + callName_ + "/" ),

    directDir_( "postProcessedWaves2Foam/" + callName_ + "/" )
{
}


postProcessingWaves::~postProcessingWaves()
{}


autoPtr<postProcessingWaves> postProcessingWaves::New
(
    const Time& rT,
    const dictionary& actionProp,
    const word& action
)
{
    postProcessingWavesConstructorTable::iterator cstrIter =
            postProcessingWavesConstructorTablePtr_->find(action);

    if (cstrIter == postProcessingWavesConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "postProcessingWaves::New(const dictionary&)"
        )   << "Unknown post processing method: " << action
            << endl << endl
            << "Valid methods are :" << endl
            << postProcessingWavesConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<postProcessingWaves>(cstrIter()( rT, actionProp, action));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
