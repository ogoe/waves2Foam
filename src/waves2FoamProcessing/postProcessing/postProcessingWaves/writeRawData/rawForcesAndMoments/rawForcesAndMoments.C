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

#include "rawForcesAndMoments.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(rawForcesAndMoments, 0);
addToRunTimeSelectionTable
(
    postProcessingWaves,
    rawForcesAndMoments,
    postProcessingWaves
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //


void rawForcesAndMoments::resizeFields
(
    List<std::pair<scalar, label> >& timeLabel,
    vectorField& field0,
    vectorField& field1,
    label N
)
{
    // Initialise timeLabel
    timeLabel.setSize(N);

    field0.setSize(N);

    field1.setSize(N);
}


void rawForcesAndMoments::writeRawData
(
    const List<std::pair<scalar, label> >& timeLabel,
    const vectorField& forces,
    const vectorField& moments
)
{
    // Write the time vector
    scalarField output0( timeLabel.size(), 0.0 );
    vectorField output1( timeLabel.size(), vector::zero );

    {
        forAll (timeLabel, labeli)
        {
            output0[labeli] = timeLabel[labeli].first;
        }

        std::stringstream ss;
        ss << callName_ << "_time";
        writeIOScalarField(output0, ss.str() );
    }

    // Write the names, indexing and dt (= -1 because of raw data format) information
    // Open a dictionary used for outputting data
    wordList wl(2);
    wl[0] = "forces";
    wl[1] = "moments";

    writeNameDict(-1, wl);

    // Write the forces as index 0
    {
        forAll (timeLabel, labeli)
        {
            output1[labeli] = forces[ timeLabel[labeli].second ];
        }

        std::stringstream ss;
        ss << callName_ << "_0";

        writeIOVectorField(output1, ss.str() );
    }

    // Write the moments as index 1
    {
        forAll (timeLabel, labeli)
        {
            output1[labeli] = moments[ timeLabel[labeli].second ];
        }

        std::stringstream ss;
        ss << callName_ << "_1";

        writeIOVectorField(output1, ss.str() );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


rawForcesAndMoments::rawForcesAndMoments
(
    const Time& rT,
    const dictionary& actionProp,
    const word& action
)
:
    postProcessingWaves( rT, actionProp, action ),

    inputDir_( actionProperties_.lookup("inputDir") ),

    removeDuplicate_( Switch(actionProperties_.lookup("removeDuplicate")) )
{
    getTimeDirs( inputDir_, timeDirs_ );
}


rawForcesAndMoments::~rawForcesAndMoments()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void rawForcesAndMoments::evaluate()
{
    List<std::pair<scalar, label> > timeLabel;
    vectorField forces, moments;

    readForceAndMomentData(timeLabel, forces, moments);

    writeRawData(timeLabel, forces, moments);
}


void rawForcesAndMoments::readForceAndMomentData
(
    List<std::pair<scalar, label> >& timeLabel,
    vectorField& forces,
    vectorField& moments
)
{
    scalar val(0.0);
    string dummy;
    label Nentries(0);

    resizeFields( timeLabel, forces, moments, 10000 );

    forAll (timeDirs_, timeI)
    {
        scalar truncateReading(0);

        if (removeDuplicate_ && timeI < timeDirs_.size() -1)
        {
            truncateReading = std::atof( timeDirs_[timeI + 1].c_str() );
        }
        else
        {
            truncateReading = GREAT;
        }

        std::stringstream ss;
        ss << inputDir_ << "/" << timeDirs_[timeI] << "/forces.dat";

        std::ifstream input;
        input.open( (ss.str()).c_str() );

        std::string line;

        // Discard the first line
        std::getline(input, line);

#if OFPLUSBRANCH==1
    // Discard two lines
    std::getline(input, line);
    std::getline(input, line);
#elif EXTBRANCH == 0
    #if OFVERSION >= 222
        // Discard yet another line
        std::getline(input, line);

        #if OFVERSION >= 230
            // Discard yet and yet another line
            std::getline(input, line);
        #endif
    #endif
#endif

        // Extracting time and overtopping flux vector
        while (std::getline( input, line ))
        {
            std::istringstream iss(line);

            // Reading the time instance
            iss >> val;

            if (truncateReading <= val)
            {
                break;
            }

            timeLabel[Nentries].first = val;
            timeLabel[Nentries].second = Nentries;

            vector temp( vector::zero );

            // Reading the first vector component with starting parenteres
            iss >> dummy;
#if EXTBRANCH==1
            temp.x() = std::atof((dummy.substr(3,dummy.size()-1)).c_str());
#elif OFPLUSBRANCH==1
            temp.x() = std::atof((dummy.substr(2,dummy.size()-1)).c_str());
#else
    #if 220<=OFVERSION
            temp.x() = std::atof((dummy.substr(3,dummy.size()-1)).c_str());
    #else
            temp.x() = std::atof((dummy.substr(2,dummy.size()-1)).c_str());
    #endif
#endif
            // Reading the second vector component. Simple scalar
            iss >> val;
            temp.y() = val;

            // Reading the third vector component with ending parenteres
            iss >> dummy;
            temp.z() = std::atof((dummy.substr(0,dummy.size()-1)).c_str());

            forces[Nentries] = temp;

            // Ignore the viscous forces
            iss >> dummy; iss >> dummy; iss >> dummy;

            // Reading the first vector component with starting parenteres
            iss >> dummy;
            temp.x() = std::atof((dummy.substr(2,dummy.size()-1)).c_str());

            // Reading the second vector component. Simple scalar
            iss >> val;
            temp.y() = val;

            // Reading the third vector component with ending parenteres
            iss >> dummy;
            temp.z() = std::atof((dummy.substr(0,dummy.size()-1)).c_str());

            moments[Nentries] = temp;

            Nentries++;

            if (Nentries == timeLabel.size())
            {
                resizeFields( timeLabel, forces, moments, 2*Nentries );
            }
        }

        input.close();
    }

    resizeFields( timeLabel, forces, moments, Nentries );

    std::sort( timeLabel.begin(), timeLabel.end(), pairSortA );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
