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

#include "ensembleAverage.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(ensembleAverage, 0);
addToRunTimeSelectionTable
(
    postProcessingWaves,
    ensembleAverage,
    postProcessingWaves
);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //


void ensembleAverage::evaluateScalar()
{
    Info << "        - Ensemble average of scalar quantities (N = dt/T = "
         << deltaT_ << "/" << period_ << " = " << N_ << ")" << endl;

    List<scalarField> input = readScalarFields( indices_ );

    List<scalarField> ensAverage( input.size() );

    scalarField count( N_ + 1, 0.0 );

    forAll (input, I)
    {
        const scalarField& field( input[I] );

        scalarField& ensAv( ensAverage[I] );
        ensAv.setSize( N_ + 1, 0.0 );

        count *= 0.0;

        forAll (field, fi)
        {
            ensAv[ fi % N_ ] += field[fi];
            count[ fi % N_ ]++;
        }

        ensAv[ N_ ] = ensAv[0];
        count[ N_ ] = count[0];

        ensAv /= count;
    }

    writeScalar( ensAverage );
}


void ensembleAverage::writeScalar
(
    const List<scalarField>& ensAv
)
{
    Info << "        - Writing ensemble average to: " << directDir_.c_str()
         << this->type() << endl;

    mkDir( directDir_ + this->type() );

    autoPtr<OFstream> spectrumPtr_;

    forAll (indices_, indexi)
    {
        std::stringstream ss;
        ss << callName_ << "_" << indices_[indexi];

        spectrumPtr_.reset
        (
            new OFstream
            (
                directDir_ + "/" + this->type() + "/" + ss.str() + "_ensAv.dat"
            )
        );

        const scalarField& data( ensAv[indexi] );

        for (int i=0; i < N_; i++)
        {
            spectrumPtr_() << static_cast<scalar>(i)*deltaT_ << tab << data[i]
                           << endl;
        }

        spectrumPtr_() << static_cast<scalar>(N_)*deltaT_ << tab
                       << data[0] << endl;
    }
}


void ensembleAverage::evaluateVector()
{
    Info << "        - Ensemble average of scalar quantities (N = dt/T = "
         << deltaT_ << "/" << period_ << " = " << N_ << ")" << endl;

    List<vectorField> input = readVectorFields( indices_ );

    List<vectorField> ensAverage( input.size() );

    scalarField count( N_ + 1, 0.0 );

    forAll (input, I)
    {
        const vectorField& field( input[I] );

        vectorField& ensAv( ensAverage[I] );
        ensAv.setSize( N_ + 1, vector::zero );

        count *= 0.0;

        forAll (field, fi)
        {
            ensAv[ fi % N_ ] += field[fi];
            count[ fi % N_ ]++;
        }

        ensAv[ N_ ] = ensAv[0];
        count[ N_ ] = count[0];

        ensAv /= count;
    }

    writeVector( ensAverage );
}


void ensembleAverage::writeVector
(
    const List<vectorField>& ensAv
)
{
    Info << "        - Writing ensemble average to: " << directDir_.c_str()
         << this->type() << endl;

    mkDir( directDir_ + this->type() );

    autoPtr<OFstream> spectrumPtr_;

    forAll (indices_, indexi)
    {
        std::stringstream ss;
        ss << callName_ << "_" << indices_[indexi];

        spectrumPtr_.reset
        (
            new OFstream
            (
                directDir_ + "/" + this->type() + "/" + ss.str() + "_ensAv.dat"
            )
        );

        const vectorField& data( ensAv[indexi] );

        for (int i=0; i < N_; i++)
        {
            spectrumPtr_() << static_cast<scalar>(i)*deltaT_
                           << tab << data[i].x()
                           << tab << data[i].y()
                           << tab << data[i].z() << endl;
        }

        spectrumPtr_() << static_cast<scalar>(N_)*deltaT_
                       << tab << data[0].x()
                       << tab << data[0].y()
                       << tab << data[0].z() << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


ensembleAverage::ensembleAverage
(
    const Time& rT,
    const dictionary& actionProp,
    const word& action
)
:
    postProcessingWaves( rT, actionProp, action ),

#   include "../../dataDict.H"
{
    readIndices( dataDict_, indices_ );

    deltaT_ = readDeltaT( dataDict_ );

    period_ = readScalar( actionProperties_.lookup("period") );

    if
    (
        Foam::mag(static_cast<label>(period_/deltaT_)*deltaT_ - period_)
        < 10*SMALL
    )
    {
        N_ = static_cast<label>( period_/deltaT_ );
    }
    else
    {
        FatalErrorIn
        (
            "ensembleAverage::ensembleAverage(const fvMesh& mesh, ...)"
        )
        << "    There is not an integer number of time steps per period.\n"
        << "    T = " << period_ << " s and " << "dt = " << deltaT_ << " s.\n"
        << "    Yielding a total of " << period_/deltaT_
        << " dt per period, which is non-integer."
        << exit(FatalError);
    }
}


ensembleAverage::~ensembleAverage()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void ensembleAverage::evaluate()
{
    if (dataType() == "scalar")
    {
        evaluateScalar();
    }
    else if (dataType() == "vector")
    {
        evaluateVector();
    }
    else
    {
        notImplemented("Only scalars and vectors are supported.");
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
