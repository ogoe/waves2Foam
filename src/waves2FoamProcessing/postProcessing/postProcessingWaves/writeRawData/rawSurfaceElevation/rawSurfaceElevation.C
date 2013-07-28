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

#include "rawSurfaceElevation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(rawSurfaceElevation, 0);
addToRunTimeSelectionTable
(
    postProcessingWaves,
    rawSurfaceElevation,
    postProcessingWaves
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //


void rawSurfaceElevation::resizeFields
(
    List<std::pair<scalar, label> >& timeLabel,
    List<scalarField>& etas,
    label N
)
{
    // Initialise timeLabel
    timeLabel.setSize(N);

    forAll (etas, etaI)
    {
        scalarField& eta( etas[etaI] );
        eta.setSize(N);
    }
}


void rawSurfaceElevation::writeRawData
(
    const List<std::pair<scalar, label> >& timeLabel,
    const scalarField& x,
    const scalarField& y,
    const scalarField& z,
    const List<scalarField>& etas
)
{
    // Write the time vector
    scalarField output( timeLabel.size(), 0.0 );

    {
        forAll (timeLabel, labeli)
        {
            output[labeli] = timeLabel[labeli].first;
        }

        std::stringstream ss;
        ss << callName_ << "_time";
        writeIOScalarField(output, ss.str() );
    }

    // Write the XYZ, indexing and dt (= -1 because of raw data format)
    // information
    writeXYZDict(-1.0, x, y, z);

    // Write the surface elevation fields
    forAll (etas, etaI)
    {
        const scalarField& eta( etas[etaI] );

        // Rearrange according to indexing
        forAll (timeLabel, labeli)
        {
            output[labeli] = eta[ timeLabel[labeli].second];
        }

        std::stringstream ss;
        ss << callName_ << "_" << etaI;

        writeIOScalarField(output, ss.str() );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


rawSurfaceElevation::rawSurfaceElevation
(
    const Time& rT,
    const dictionary& actionProp,
    const word& action
)
:
    postProcessingWaves( rT, actionProp, action ),

    inputDir_
    (
        actionProperties_.lookupOrDefault<word>("inputDir", "surfaceElevation")
    ),

    removeDuplicate_( Switch(actionProperties_.lookup("removeDuplicate")) )
{
    getTimeDirs(inputDir_, timeDirs_);
}


rawSurfaceElevation::~rawSurfaceElevation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void rawSurfaceElevation::evaluate()
{
    List<std::pair<scalar, label> > timeLabel;
    scalarField x, y, z;
    List<scalarField> etas;

    readSurfaceElevationData(timeLabel, x, y, z, etas);

    writeRawData(timeLabel, x, y, z, etas);
}


void rawSurfaceElevation::readSurfaceElevationData
(
    List<std::pair<scalar, label> >& timeLabel,
    scalarField& x,
    scalarField& y,
    scalarField& z,
    List<scalarField>& etas
)
{
    label Ngauges(0);
    scalar val(0.0);
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
        ss << inputDir_ << "/" << timeDirs_[timeI] << "/surfaceElevation.dat";

        std::ifstream input;
        input.open( (ss.str()).c_str() );

        std::string line;

        if (timeI == 0)
        {
            // Reading the x-coordinates
            {
                // Discard the first line
                std::getline( input, line);
                std::getline( input, line);

                std::istringstream iss(line);

                // Discard the first data entry
                iss >> val;

                while (iss >> val)
                {
                    x.setSize( Ngauges + 1 );
                    x[Ngauges++] = val;
                }

                y.setSize( Ngauges );
                z.setSize( Ngauges );
            }

            // Reading the y-coordinates
            {
                std::getline( input, line);

                std::istringstream iss(line);

                iss >> val;

                for (int i=0; i < Ngauges; i++)
                {
                    iss >> val;
                    y[i] = val;
                }
            }

            // Reading the z-coordinates
            {
                std::getline( input, line);

                std::istringstream iss(line);

                iss >> val;

                for (int i=0; i < Ngauges; i++)
                {
                    iss >> val;
                    z[i] = val;
                }
            }

            etas.setSize( Ngauges );

            resizeFields( timeLabel, etas, 10000 );
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
            iss >> val;

            if (truncateReading <= val)
            {
                break;
            }

            timeLabel[Nentries].first = val;
            timeLabel[Nentries].second = Nentries;

            forAll (etas, etaI)
            {
                iss >> val;
                scalarField& eta( etas[ etaI ] );

                eta[Nentries] = val;
            }

            Nentries++;

            if (Nentries == timeLabel.size())
            {
                resizeFields( timeLabel, etas, 2*Nentries );
            }
        }

        input.close();
    }

    resizeFields( timeLabel, etas, Nentries );

    std::sort( timeLabel.begin(), timeLabel.end(), pairSortA );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
