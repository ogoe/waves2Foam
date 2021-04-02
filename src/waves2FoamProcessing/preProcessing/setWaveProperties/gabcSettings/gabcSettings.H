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

Class
    Foam::gabcSettings

Description


SourceFiles
    gabcSettings.C

Author
    Niels Gjoel Jacobsen, Deltares.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef gabcSettings_H
#define gabcSettings_H

#include "fvMesh.H"
#include "fvCFD.H"
#include "volFields.H"
#include "IOdictionary.H"

#include "OFstream.H"

#include "autoPtr.H"
#include "runTimeSelectionTables.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class gabcSettings Declaration
\*---------------------------------------------------------------------------*/

class gabcSettings
{

protected:

    // Protected data
    dictionary& dict_;
    scalar depth_;

    label Nspaces_;

    void writeGiven( Ostream& , word );

    void writeDerived( Ostream&, word, scalar);

    void writeDerived( Ostream&, word, word);

    void addITstream(Ostream&, const word&, const ITstream&);

    inline void writeDepth(Ostream& os)
    {
        if (depth_ < 0)
        {
            writeGiven(os, "depth");
        }
        else
        {
            writeDerived(os, "depth", depth_);
        }
    }

private:

    // Private Member Functions
        void lineFormatting( Ostream&, const word&);

        //- Disallow default bitwise copy construct
        gabcSettings(const gabcSettings&);

        //- Disallow default bitwise assignment
        void operator=(const gabcSettings&);

public:

    //- Runtime type information
    TypeName("gabcSettings");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            gabcSettings,
            gabcSettings,
            (
                dictionary& dict, scalar depth
            ),
            (dict, depth)
        );

    // Constructors

        //- Construct from components
        gabcSettings
        (
            dictionary&,
            scalar
        );

    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<gabcSettings> New
        (
            dictionary&,
            scalar
        );

    // Destructor

        virtual ~gabcSettings();

    // Member Functions

        virtual void set(Ostream&) = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
