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

#include "rawOvertopping.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(rawOvertopping, 0);
addToRunTimeSelectionTable
(
    postProcessingWaves,
    rawOvertopping,
    postProcessingWaves
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //


void rawOvertopping::resizeFields
(
    List<std::pair<scalar, label> >& timeLabel,
    List<vectorField>& Us,
    label N
)
{
    // Initialise timeLabel
    timeLabel.setSize(N);

    forAll (Us, UI)
    {
        vectorField& U( Us[UI] );
        U.setSize(N);
    }
}


void rawOvertopping::writeRawData
(
    const List<std::pair<scalar, label> >& timeLabel,
    const wordList& OTnames,
    const List<vectorField>& OTs
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

    // Write the names, indexing and dt (= -1 because of raw data format)
    // information. Open a dictionary used for outputting data
    writeNameDict(-1, OTnames);

    // Write the surface elevation fields
    forAll (OTs, OTI)
    {
        const vectorField& OT( OTs[OTI] );

        // Rearrange according to indexing
        forAll (timeLabel, labeli)
        {
            output1[labeli] = OT[ timeLabel[labeli].second ];
        }

        std::stringstream ss;
        ss << callName_ << "_" << OTI;

        writeIOVectorField(output1, ss.str() );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


rawOvertopping::rawOvertopping
(
    const Time& rT,
    const dictionary& actionProp,
    const word& action
)
:
    postProcessingWaves( rT, actionProp, action ),

    inputDir_
    (
        actionProperties_.lookupOrDefault<word>("inputDir", "overtopping")
    ),

    scaleFlux_
    (
        actionProperties_.lookupOrDefault<scalar>("overtoppingScaling",1.0)
    ),

    removeDuplicate_( Switch(actionProperties_.lookup("removeDuplicate")) )
{
    getTimeDirs(inputDir_, timeDirs_);
}


rawOvertopping::~rawOvertopping()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void rawOvertopping::evaluate()
{
    List<std::pair<scalar, label> > timeLabel;
    wordList OTnames;
    List<vectorField> OTs;

    readOvertoppingData(timeLabel, OTnames, OTs);

    writeRawData(timeLabel, OTnames, OTs);
}


void rawOvertopping::readOvertoppingData
(
    List<std::pair<scalar, label> >& timeLabel,
    wordList& OTnames,
    List<vectorField>& OTs
)
{
    label Nprobes(0);
    scalar val(0.0);
    string dummy;
    label Nentries(0);

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
        ss << inputDir_ << "/" << timeDirs_[timeI] << "/overtopping.dat";

        std::ifstream input;
        input.open( (ss.str()).c_str() );

        std::string line;

        if (timeI == 0)
        {
            std::getline( input, line);

            std::istringstream iss(line);

            // Discard first string
            iss >> dummy;

            while (iss >> dummy)
            {
                OTnames.setSize( Nprobes + 1 );
                OTnames[Nprobes++] = dummy;
            }

            OTs.setSize( Nprobes );

            resizeFields( timeLabel, OTs, 10000 );
        }
        else
        {
            // Disgard the first line
            std::getline( input, line);
        }

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

            forAll (OTs, OTI)
            {
                vector temp( vector::zero );

                // Reading the first vector component with starting parenteres
                iss >> dummy;
                temp.x() = std::atof((dummy.substr(1,dummy.size()-1)).c_str());

                // Reading the second vector component. Simple scalar
                iss >> val;
                temp.y() = val;

                // Reading the third vector component with ending parenteres
                iss >> dummy;
                temp.z() = std::atof((dummy.substr(0,dummy.size()-1)).c_str());

                vectorField& OT( OTs[ OTI ] );

                OT[Nentries] = scaleFlux_*temp;
            }

            Nentries++;

            if (Nentries == timeLabel.size())
            {
                resizeFields( timeLabel, OTs, 2*Nentries );
            }
        }

        input.close();
    }

    resizeFields( timeLabel, OTs, Nentries );

    std::sort( timeLabel.begin(), timeLabel.end(), pairSortA );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
