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

#include "rawVelocityProbes.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(rawVelocityProbes, 0);
addToRunTimeSelectionTable
(
    postProcessingWaves,
    rawVelocityProbes,
    postProcessingWaves
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //


void rawVelocityProbes::resizeFields
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


void rawVelocityProbes::writeRawData
(
    const List<std::pair<scalar, label> >& timeLabel,
    const scalarField& x,
    const scalarField& y,
    const scalarField& z,
    const List<vectorField>& Us
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

    // Write the XYZ, indexing and dt (= -1 because of raw data format)
    // information
    writeXYZDict(-1.0, x, y, z);

    // Write the surface elevation fields
    forAll (Us, UI)
    {
        const vectorField& U( Us[UI] );

        // Rearrange according to indexing
        forAll (timeLabel, labeli)
        {
            output1[labeli] = U[ timeLabel[labeli].second ];
        }

        std::stringstream ss;
        ss << callName_ << "_" << UI;

        writeIOVectorField(output1, ss.str() );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


rawVelocityProbes::rawVelocityProbes
(
    const Time& rT,
    const dictionary& actionProp,
    const word& action
)
:
    postProcessingWaves( rT, actionProp, action ),

    inputDir_( actionProperties_.lookupOrDefault<word>("inputDir", "probes") ),

    removeDuplicate_( Switch(actionProperties_.lookup("removeDuplicate")) ),

    R_(1,0,0,0,1,0,0,0,1)
{
    getTimeDirs(inputDir_, timeDirs_);

    if (actionProperties_.lookupOrDefault<Switch>("rotateVelocity",false))
    {
#if OFVERSION<230
        coordinateRotation cr( actionProperties_.subDict("rotate") );

        R_ = cr.R().T();
#else
        autoPtr<coordinateRotation> cr
            (
                Foam::coordinateRotation::New
                (
                    actionProperties_.subDict("rotate")
                )
            );

        R_ = cr->R().T();
#endif
    }

    Info << "       - Uses the rotation vector: " << R_ << endl;
}


rawVelocityProbes::~rawVelocityProbes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void rawVelocityProbes::evaluate()
{
    List<std::pair<scalar, label> > timeLabel;
    scalarField x, y, z;
    List<vectorField> Us;

    readVelocityProbeData(timeLabel, x, y, z, Us);

    writeRawData(timeLabel, x, y, z, Us);
}


void rawVelocityProbes::readVelocityProbeData
(
    List<std::pair<scalar, label> >& timeLabel,
    scalarField& x,
    scalarField& y,
    scalarField& z,
    List<vectorField>& Us
)
{
    label Nprobes(0);
    scalar val(0.0);
    string dummy;
    label Nentries(0);

    forAll (timeDirs_, timeI)
    {
        scalar truncateReading(0);

        if (removeDuplicate_ && timeI < timeDirs_.size() - 1)
        {
            truncateReading = std::atof( timeDirs_[timeI + 1].c_str() );
        }
        else
        {
            truncateReading = GREAT;
        }

        std::stringstream ss;
        ss << inputDir_ << "/" << timeDirs_[timeI] << "/U";

        std::ifstream input;
        input.open( (ss.str()).c_str() );

        std::string line;

        if (timeI == 0)
        {
            // Reading the x-coordinates
            {
                std::getline( input, line);

                std::istringstream iss(line);

                // Discard the first data entry
                iss >> dummy;
                iss >> dummy;

                while (iss >> val)
                {
                    x.setSize( Nprobes + 1 );
                    x[Nprobes++] = val;
                }

                y.setSize( Nprobes );
                z.setSize( Nprobes );
            }

            // Reading the y-coordinates
            {
                std::getline( input, line);

                std::istringstream iss(line);

                iss >> dummy; iss >> dummy;

                for (int i=0; i < Nprobes; i++)
                {
                    iss >> val;
                    y[i] = val;
                }
            }

            // Reading the z-coordinates
            {
                std::getline( input, line);

                std::istringstream iss(line);

                iss >> dummy; iss >> dummy;

                for (int i=0; i < Nprobes; i++)
                {
                    iss >> val;
                    z[i] = val;
                }
            }

            // Discard the fourth line
            std::getline( input, line);

            Us.setSize( Nprobes );

            resizeFields( timeLabel, Us, 10000 );
        }
        else
        {
            // Disgard the first four lines
            for (int i=0; i < 4; i++)
            {
                std::getline( input, line);
            }
        }

        // Extracting time and surface elevation
        while (std::getline( input, line ))
        {
            std::istringstream iss(line);

            // Reading the time instance
            iss >> val;

            if (truncateReading <= val)
                break;

            timeLabel[Nentries].first = val;
            timeLabel[Nentries].second = Nentries;

            forAll (Us, UI)
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

                vectorField& U( Us[ UI ] );

                // Rotates the velocities into a given local coordinate system
                U[Nentries] = (R_ & temp);
            }

            Nentries++;

            if (Nentries == timeLabel.size())
            {
                resizeFields( timeLabel, Us, 2*Nentries );
            }
        }

        input.close();
    }

    resizeFields( timeLabel, Us, Nentries );

    std::sort( timeLabel.begin(), timeLabel.end(), pairSortA );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
