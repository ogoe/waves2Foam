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

#include "zeroCrossing.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(zeroCrossing, 0);
addToRunTimeSelectionTable
(
    postProcessingWaves,
    zeroCrossing,
    postProcessingWaves
);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //


void zeroCrossing::evaluateScalar()
{
    Info << "        - Zero-crossing analysis. Crossing method is: "
         << crossingType_ << endl;

    // Get the time axis
    std::stringstream ss;
    ss << callName_ << "_time";

    scalarField t;
    t = readIOScalarField(ss.str());

    // Get the surface elevation fields
    List<scalarField> input = readScalarFields(indices_);

    // Prepare the output data
    List<scalarField> startingTime(input.size());
    List<scalarField> duration(input.size());
    List<scalarField> maxValue(input.size());
    List<scalarField> minValue(input.size());

    if (crossingType_ == "up")
    {
        crossingAnalysis
        (
            t,
            input,
            startingTime,
            duration,
            maxValue,
            minValue,
            2
        );

        writeScalar(startingTime, duration, maxValue, minValue, crossingType_);
    }
    else if (crossingType_ == "down")
    {
        crossingAnalysis
        (
            t,
            input,
            startingTime,
            duration,
            maxValue,
            minValue,
            -2
        );

        writeScalar(startingTime, duration, maxValue, minValue, crossingType_);
    }
    else if (crossingType_ == "both")
    {
        // First the up-crossing analysis
        crossingAnalysis
        (
            t,
            input,
            startingTime,
            duration,
            maxValue,
            minValue,
            2
        );

        writeScalar(startingTime, duration, maxValue, minValue, "up");

        // Then the down-crossing analysis
        crossingAnalysis
        (
            t,
            input,
            startingTime,
            duration,
            maxValue,
            minValue,
            -2
        );

        writeScalar(startingTime, duration, maxValue, minValue, "down");
    }
    else
    {
        FatalErrorIn("void zeroCrossing::evaluateScalar()")
            << "\nDefine either of the following crossing methods for the\n"
            << "zero-crossing analysis: up, down or both\n"
            << endl << endl << exit(FatalError);
    }
}


void zeroCrossing::writeScalar
(
    const List<scalarField>& startingTime,
    const List<scalarField>& duration,
    const List<scalarField>& maxValue,
    const List<scalarField>& minValue,
    const word lastName
) const
{
    Info << "        - Writing zero-crossing results to: " << directDir_.c_str()
         << this->type() << endl;

    mkDir(directDir_ + this->type());

    autoPtr<OFstream> crossingPtr_;

    forAll (indices_, indexi)
    {
        std::stringstream ss;
        ss << callName_ << "_" << lastName << "_" << indices_[indexi];

        crossingPtr_.reset
        (
            new OFstream
            (
                directDir_ + "/" + this->type() + "/" + ss.str() + ".dat"
            )
        );

        const scalarField& sT = startingTime[indexi];
        const scalarField& d = duration[indexi];
        const scalarField& xV = maxValue[indexi];
        const scalarField& nV = minValue[indexi];

        for (int i = 0; i < sT.size(); i++)
        {
            crossingPtr_() << sT[i] << tab << d[i] << tab
                           << nV[i] << tab << xV[i] << endl;
        }
    }
}


void zeroCrossing::crossingAnalysis
(
    const scalarField& t,
    const List<scalarField>& input,
    List<scalarField>& startingTime,
    List<scalarField>& duration,
    List<scalarField>& maxValue,
    List<scalarField>& minValue,
    const label crossingVal
) const
{
    forAll (input, indexi)
    {
        // Initialise the output fields
        scalarField& sT = startingTime[indexi];
        sT.setSize(500, 0);

        scalarField& d = duration[indexi];
        d.setSize(500, 0);

        scalarField& xV = maxValue[indexi];
        xV.setSize(500, 0);

        scalarField& nV = minValue[indexi];
        nV.setSize(500, 0);

        // Get reference to input field and create crossing field
        const scalarField& inp = input[indexi];
        scalarField searchField = Foam::pos(inp) - Foam::neg(inp);

        label count = -1;

        scalar MAX = -GREAT;
        scalar MIN = GREAT;
        scalar tStart = t[0];

        for (label n = 0; n < searchField.size() - 1; n++)
        {
            if (inp[n] > MAX)
            {
                MAX = inp[n];
            }

            if (inp[n] < MIN)
            {
                MIN = inp[n];
            }

            if (Foam::mag(searchField[n + 1] - searchField[n] - crossingVal) < SMALL)
            {
                scalar tNew = t[n] - inp[n]*deltaT_/(inp[n + 1] - inp[n]);

                if (count > -1)
                {
                    sT[count] = tStart;
                    d[count] = tNew - tStart;
                    xV[count] = MAX;
                    nV[count] = MIN;
                }

                MAX = -GREAT;
                MIN = GREAT;
                tStart = tNew;

                count++;

                if (sT.size() == count)
                {
                    sT.setSize(2*count);
                    d.setSize(2*count);
                    xV.setSize(2*count);
                    nV.setSize(2*count);
                }
            }
        }

        sT.setSize(count);
        d.setSize(count);
        xV.setSize(count);
        nV.setSize(count);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


zeroCrossing::zeroCrossing
(
    const Time& rT,
    const dictionary& actionProp,
    const word& action
)
:
    postProcessingWaves( rT, actionProp, action ),

    crossingType_(actionProp.lookup("crossingType")),

#   include "dataDict.H"
{
    readIndices(dataDict_, indices_);

    deltaT_ = readDeltaT(dataDict_);
}


zeroCrossing::~zeroCrossing()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void zeroCrossing::evaluate()
{
    if (dataType() == "scalar")
    {
        evaluateScalar();
    }
    else
    {
        notImplemented("Only scalars are supported.");
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
