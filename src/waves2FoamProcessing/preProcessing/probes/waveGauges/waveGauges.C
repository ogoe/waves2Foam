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

#include "waveGauges.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(waveGauges, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void waveGauges::writeVTKFormat
(
    const word& name,
    const pointField& pp,
    const point& addPoint
)
{
    autoPtr<OFstream> vtk;

    // Writing the lines as VTK format
    vtk.reset(new OFstream("waveGaugesNProbes/" + name + ".vtk" ));

    // Writing header
    vtk() << "# vtk DataFile Version 3.0" << nl << "vtk output" << nl
          << "ASCII" << nl << "DATASET POLYDATA" << endl;

    // Writing points
    vtk() << "POINTS " << 2*pp.size() << " float" << endl;

    forAll(pp, pointi)
    {
        point p(pp[pointi]);
        vtk() << p.x() << " " << p.y() << " " << p.z() << endl;
    }

    forAll(pp, pointi)
    {
        point p(pp[pointi] + addPoint);
        vtk() << p.x() << " " << p.y() << " " << p.z() << endl;
    }

    // Writing lines
    vtk() << "LINES " << pp.size() << " " << 3*pp.size() << endl;
    forAll (pp, pointi)
    {
        vtk() << "2 " << pointi << " " << pointi + pp.size() << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


waveGauges::waveGauges
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),

    gaugeDict_(dict)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void waveGauges::evaluate(const word& name)
{
    point addPoint(gaugeDict_.lookup("add"));
    word  vertAxis(word(gaugeDict_.lookup("axis")));

    autoPtr<Foam::pointDistributions> pd
    (
        Foam::pointDistributions::New(mesh_, gaugeDict_)
    );

    pointField pp(pd->evaluate());

    autoPtr<OFstream> gauges;

    // Writing the sets file
    gauges.reset(new OFstream("waveGaugesNProbes/" + name + "_sets"));

    gauges() << "sets" << nl << token::BEGIN_LIST << nl << incrIndent;

    forAll (pp, pointi)
    {
        gauges() << indent << "gauge_" << pointi << nl << indent
                 << token::BEGIN_BLOCK << incrIndent << nl;
        gauges() << indent << "type         face"
                 << token::END_STATEMENT << nl;
        gauges() << indent << "axis         " << vertAxis
                 << token::END_STATEMENT << nl;
        gauges() << indent << "start        " << pp[pointi]
                 << token::END_STATEMENT << nl;
        gauges() << indent << "end          " << pp[pointi] + addPoint
                 << token::END_STATEMENT << nl;
        gauges() << indent << "nPoints      100" << token::END_STATEMENT << nl;
        gauges() << decrIndent << indent << token::END_BLOCK << nl << nl;
    }

    gauges() << decrIndent << token::END_LIST << token::END_STATEMENT << nl;

    // Writing the file to be included in controlDict
    gauges.reset(new OFstream("waveGaugesNProbes/" + name + "_controlDict"));

    gauges() << incrIndent << indent << name << nl << indent
             << token::BEGIN_BLOCK << nl << incrIndent;
    gauges() << indent << "type               surfaceElevation;" << nl;
    gauges() << indent << "functionObjectLibs ( \"libwaves2Foam.so\" );" << nl;
    gauges() << nl;
#if OFPLUSBRANCH==1 && 1606<OFVERSION
    word wc(gaugeDict_.lookup("writeControl"));
    scalar wi = readScalar(gaugeDict_.lookup("writeInterval"));
    gauges() << indent << "writeControl       " << wc << ";" << nl;
    gauges() << indent << "writeInterval      " << wi << ";" << nl;
#else
    gauges() << indent << "outputControl      timeStep;"
             << " // Alternative: outputTime" << nl;
    gauges() << indent << "outputInterval      1;" << nl << nl;
    gauges() << indent << "//Additional output controls in waves2Foam" << nl;
    gauges() << indent << "//samplingStartTime  -1;" << nl;
    gauges() << indent << "//surfaceSampleDeltaT 0.025;" << nl;
#endif
    gauges() << nl;
    gauges() << indent << "setFormat          raw;" << nl;
    gauges() << indent << "interpolationScheme cellPointFace;" << nl;
    gauges() << indent << "fields (" << Foam::waves2Foam::aName() << ");" 
             << nl << nl;
    gauges() << indent << "#includeIfPresent \"../waveGaugesNProbes/" << name
             << "_sets\";" << nl << nl;
    gauges() << decrIndent << indent << token::END_BLOCK << decrIndent << nl;

    // Writing the surfaceElevationDict
    gauges.reset
    (
        new OFstream("waveGaugesNProbes/" + name + "surfaceElevationDict")
    );

    mesh_.writeBanner(gauges());

    // Write the file information. Class name is not correct when
    // using wOut.writeHeader( os ); hence manual entries
    gauges() << "FoamFile" << nl;
    gauges() << token::BEGIN_BLOCK << incrIndent << nl;
    gauges() << indent << "version" << tab << IOstream::currentVersion
             << token::END_STATEMENT << nl;
    gauges() << indent << "format" << tab << "ascii;" << nl;
    gauges() << indent << "class" << tab << "dictionary;" << nl;
    gauges() << indent << "object" << tab << "surfaceElevationDict;" << nl;
    gauges() << decrIndent << indent << token::END_BLOCK << nl;

    // Write the divider
    mesh_.writeDivider(gauges());
    gauges() << nl << nl;

    gauges() << "setFormat           raw;" << nl;
    gauges() << "interpolationScheme cellPointFace;" << nl;
    gauges() << indent << "fields (" << Foam::waves2Foam::aName() << ");" 
             << nl << nl;
    gauges() << "#includeIfPresent  \"../waveGaugesNProbes/" << name
             << "_sets\";" << nl;

    mesh_.writeEndDivider(gauges());

    if (gaugeDict_.lookupOrDefault<Switch>("writeVTK", true))
    {
        writeVTKFormat(name, pp, addPoint);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
