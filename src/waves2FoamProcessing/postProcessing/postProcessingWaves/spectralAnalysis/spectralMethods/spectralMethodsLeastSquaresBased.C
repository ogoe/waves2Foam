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

#include "spectralMethodsLeastSquaresBased.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(spectralMethodsLeastSquaresBased, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void spectralMethodsLeastSquaresBased::computePowerSpectrum
(
    const scalarField& t,
    const scalarField& input,
    const label& N,
    const scalar& f,
    scalarField& bOut
)
{
    // A cosine and a sine for each frequency and
    // one output for the mean value
    bOut.setSize( input.size() );
    bOut = input;

    List<scalarField> matrix( 2*N + 1 );

    // Set ones in the right most column
    scalarField& m( matrix[2*N] );
    m.setSize( t.size(), 1.0 );

    // Set the cosine components - every second
    // column starting from column index 0
    for (int i = 0; i < N; i++)
    {
        scalarField& m( matrix[ 2*i] );
        m.setSize( t.size(), 0.0);

        m = Foam::cos( static_cast<scalar>(i + 1)*2 * M_PI*f * t );
    }

    // Set the sine components - every second
    // column starting from column index 0
    for (int i = 0; i < N; i++)
    {
        scalarField& m( matrix[ 2*i + 1] );
        m.setSize( t.size(), 0.0);

        m = Foam::sin( static_cast<scalar>(i + 1)*2 * M_PI*f * t );
    }

    solve( matrix, bOut );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


spectralMethodsLeastSquaresBased::spectralMethodsLeastSquaresBased
(
    const Time& rT,
    const dictionary& actionProp
)
{
}


spectralMethodsLeastSquaresBased::~spectralMethodsLeastSquaresBased()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void spectralMethodsLeastSquaresBased::solve
(
    const List<scalarField>& A,
    scalarField& b
)
{
    // Create the least-squares right and left hand sides
    label N( A.size() );

#if EXTBRANCH==1
    scalarSquareMatrix AtA(N, 0.0);
#elif OFPLUSBRANCH==1
    SquareMatrix<scalar> AtA(N, N);
#else
    SquareMatrix<scalar> AtA(N, N);
#endif
    scalarField Atb( N, 0.0 );

    // Fill the matrix elements
    for (int i=0; i<N; i++)
    {
        const scalarField& ai( A[i] );

        for (int j=0; j<N; j++)
        {
            const scalarField& aj( A[j] );
            AtA[i][j] = Foam::sum( ai*aj );
        }

        Atb[i] = Foam::sum(ai*b);
    }

    // Solve the square system
#if EXTBRANCH==1
    Foam::scalarSquareMatrix::LUsolve(AtA, Atb);
#elif OFPLUSBRANCH==1
    Foam::LUsolve(AtA, Atb);
#else
    Foam::LUsolve(AtA, Atb);
#endif

    // Return solution in input field
    b = Atb;
}


scalarField spectralMethodsLeastSquaresBased::frequencies
(
    const label& N
)
{
    scalarField res(2*N + 1, 0);

    for (int i=0; i<N; i++)
    {
        res[ 2*i    ] = i + 1;
        res[ 2*i + 1] = i + 1;
    }

    return res;
}


scalarField spectralMethodsLeastSquaresBased::powerSpectra
(
    const scalarField& t,
    const scalarField& input,
    const label& N,
    const scalar& f
)
{
    scalarField res;

    computePowerSpectrum( t, input, N, f, res );

    return res;
}


List<scalarField> spectralMethodsLeastSquaresBased::powerSpectra
(
    const scalarField& t,
    const List<scalarField>& input,
    const label& N,
    const scalar& f
)
{
    List<scalarField> res( input.size() );

    forAll (input, inputi)
    {
        const scalarField i( input[inputi] );
        scalarField& r( res[inputi] );

        r = powerSpectra(t, i, N, f );
    }

    return res;
}


List<vectorField> spectralMethodsLeastSquaresBased::powerSpectra
(
    const scalarField& t,
    const List<vectorField>& input,
    const label& N,
    const scalar& f
)
{
    List<vectorField> res( input.size() );

    forAll (input, inputi)
    {
        const vectorField i( input[inputi] );
        vectorField& r( res[inputi] );
        r.setSize( 2*N + 1);

        for (int j = 0; j < 3; j++)
        {
            r.replace( j, powerSpectra(t, i.component(j), N, f ));
        }
    }

    return res;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
