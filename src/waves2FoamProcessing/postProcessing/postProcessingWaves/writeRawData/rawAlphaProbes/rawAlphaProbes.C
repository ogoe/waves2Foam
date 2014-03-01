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

#include "rawAlphaProbes.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(rawAlphaProbes, 0);
addToRunTimeSelectionTable
(
    postProcessingWaves,
    rawAlphaProbes,
    postProcessingWaves
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //


void rawAlphaProbes::resizeFields
(
    List<std::pair<scalar, label> >& timeLabel,
    List<scalarField>& alphas,
    label N
)
{
    // Initialise timeLabel
    timeLabel.setSize(N);

    forAll (alphas, alphaI)
    {
        scalarField& alpha( alphas[alphaI] );
        alpha.setSize(N);
    }
}


void rawAlphaProbes::writeRawData
(
    const List<std::pair<scalar, label> >& timeLabel,
    const scalarField& x,
    const scalarField& y,
    const scalarField& z,
    const List<scalarField>& alphas
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
    forAll (alphas, alphaI)
    {
        const scalarField& alpha( alphas[alphaI] );

        // Rearrange according to indexing
        forAll (timeLabel, labeli)
        {
            output[labeli] = alpha[ timeLabel[labeli].second];
        }

        std::stringstream ss;
        ss << callName_ << "_" << alphaI;

        writeIOScalarField(output, ss.str() );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


rawAlphaProbes::rawAlphaProbes
(
    const Time& rT,
    const dictionary& actionProp,
    const word& action
)
:
    postProcessingWaves( rT, actionProp, action ),

    inputDir_( actionProperties_.lookupOrDefault<word>("inputDir", "probes") ),

    removeDuplicate_( Switch(actionProperties_.lookup("removeDuplicate")) )
{
    getTimeDirs(inputDir_, timeDirs_);
}


rawAlphaProbes::~rawAlphaProbes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void rawAlphaProbes::evaluate()
{
    List<std::pair<scalar, label> > timeLabel;
    scalarField x, y, z;
    List<scalarField> alphas;

    readAlphaProbesData(timeLabel, x, y, z, alphas);

    writeRawData(timeLabel, x, y, z, alphas);
}


void rawAlphaProbes::readAlphaProbesData
(
    List<std::pair<scalar, label> >& timeLabel,
    scalarField& x,
    scalarField& y,
    scalarField& z,
    List<scalarField>& alphas
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
        ss << inputDir_ << "/" << timeDirs_[timeI] << "/" 
           << Foam::waves2Foam::aName();

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

            alphas.setSize( Nprobes );

            resizeFields( timeLabel, alphas, 10000 );
        }
        else
        {
            // Disgard the first four lines
            for (int i=0; i < 4; i++)
            {
                std::getline( input, line);
            }
        }

        // Extracting time and alpha
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

            forAll (alphas, alphaI)
            {
                iss >> val;
                scalarField& alpha( alphas[ alphaI ] );

                alpha[Nentries] = val;
            }

            Nentries++;

            if (Nentries == timeLabel.size())
            {
                resizeFields( timeLabel, alphas, 2*Nentries );
            }
        }

        input.close();
    }

    resizeFields( timeLabel, alphas, Nentries );

    std::sort( timeLabel.begin(), timeLabel.end(), pairSortA );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
