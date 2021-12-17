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

#include "fixedDischarge.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(fixedDischargeTrial, 0);
addToRunTimeSelectionTable(waveTheory, fixedDischargeTrial, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


fixedDischargeTrial::fixedDischargeTrial
(
    const word& subDictName,
    const fvMesh& mesh_
)
:
    waveTheory(subDictName, mesh_),

    fvMesh_(mesh_),

    boundaryName_(subDictName),

    Tsoft_(readScalar(coeffDict_.lookup("Tsoft"))),

    Qm3s_
    (
        readScalar(coeffDict_.lookup("Qm3s"))
    ),

    bottomLevel_
    (
        readScalar(coeffDict_.lookup("bottomLevel"))
    ),

    width_
    (
        readScalar(coeffDict_.lookup("domainWidth"))
    ),


    timeIndex_(-1),

    localSeaLevel_(seaLevel_),

    filterTime_(readScalar(coeffDict_.lookup("filterTime"))),

    filterDt_(filterTime_/50.0),

    filterLastSample_(-1),

    surfaceTime_(50, mesh_.time().startTime().value()),

    surfaceElevation_(50, seaLevel_)
{
}


void fixedDischargeTrial::printCoeffs()
{
    Info << "Loading wave theory: " << typeName << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fixedDischargeTrial::updateLocalSeaLevel() const
{
    if (timeIndex_ < fvMesh_.time().timeIndex())
    {
        timeIndex_ = fvMesh_.time().timeIndex();

        label patchID = fvMesh_.boundaryMesh().findPatchID(boundaryName_);
        const volScalarField& alpha
        (
            fvMesh_.thisDb().lookupObject<volScalarField>("alpha.water")
        );

        const vectorField& Sf = fvMesh_.Sf().boundaryField()[patchID];
        scalar waterArea = Foam::gSum(Foam::mag(Sf)*alpha.boundaryField()[patchID]);
        
        localSeaLevel_ = waterArea/width_ + bottomLevel_;

        if (0 < filterTime_)
        {
            // Update filtering values
            label N = surfaceTime_.size();
            if (filterLastSample_ + filterDt_ < fvMesh_.time().time().value())
            {
                Info << "IN LOOP: " << filterLastSample_ << endl;

                filterLastSample_ = fvMesh_.time().time().value();

                // Move all values one back in the arracy
                for (long index = 1; index < N; index++)
                {
                    surfaceTime_[N - index] = surfaceTime_[N - index - 1];
                    surfaceElevation_[N - index] = surfaceElevation_[N - index - 1];
                }

                surfaceTime_[0] = filterLastSample_;
                surfaceElevation_[0] = localSeaLevel_;
            }

            // Calculate the average surface elevation
            scalar averageSurfaceElevation(0);

            averageSurfaceElevation += (fvMesh_.time().time().value() - surfaceTime_[0])
                    *(localSeaLevel_ + surfaceElevation_[0]);

            for (long index = 1; index < N; index++)
            {
                averageSurfaceElevation += (surfaceTime_[index - 1] - surfaceTime_[index])
                        *(surfaceElevation_[index] + surfaceElevation_[index - 1]);
            }

            averageSurfaceElevation /= (2*(fvMesh_.time().time().value() - surfaceTime_[surfaceTime_.size() - 1]));

            // Set the local sea level to the averaged value
            localSeaLevel_ = averageSurfaceElevation;
        }
    }
}



scalar fixedDischargeTrial::factor(const scalar& time) const
{
    scalar factor(1.0);
    if (Tsoft_ > 0.0)
    {
        factor = Foam::sin(2*PI_/(4.0*Tsoft_)*Foam::min(Tsoft_, time));
    }

    return factor;
}


scalar fixedDischargeTrial::eta
(
    const point& x,
    const scalar& time
) const
{
    updateLocalSeaLevel();

    scalar eta = localSeaLevel_;

    return eta;
}


scalar fixedDischargeTrial::pExcess
(
    const point& x,
    const scalar& time
) const
{
    updateLocalSeaLevel();

    scalar result = rhoWater_*Foam::mag(g_)*localSeaLevel_;

    result += rhoWater_*(referenceLevel_ & g_);

    return result;
}


vector fixedDischargeTrial::U
(
    const point& x,
    const scalar& time
) const
{
    // Get a reference to the velocity field
    const volScalarField& alpha
        (
            fvMesh_.thisDb().lookupObject<volScalarField>(Foam::waves2Foam::aName())
        );

    // Get a reference to the patch ID and the face centres
    label patchID = fvMesh_.boundaryMesh().findPatchID(boundaryName_);
    const vectorField& Sf = fvMesh_.Sf().boundaryField()[patchID];

    // Alpha area
    scalar alphaArea = Foam::sum(alpha.boundaryField()[patchID]*Foam::mag(Sf));

    // Flux is normal to the boundary
    vector nf = -Sf[0]/Foam::mag(Sf[0]);

    // Set to 0 in case of no wetting (effectively encountered under some
    // cases with initiation).
    if (alphaArea < SMALL)
    {
        return 0*nf;
    }
    else
    {
    	return Qm3s_/alphaArea*nf*factor(time);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
