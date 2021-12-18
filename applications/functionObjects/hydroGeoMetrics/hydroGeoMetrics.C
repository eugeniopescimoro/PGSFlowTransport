/*---------------------------------------------------------------------------*\
=========                 |
\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
\\    /   O peration     | Website:  https://openfoam.org
\\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
\\/     M anipulation  |
-------------------------------------------------------------------------------
License
This file is part of OpenFOAM.

OpenFOAM is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "hydroGeoMetrics.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace functionObjects
    {
        defineTypeNameAndDebug(hydroGeoMetrics, 0);
        addToRunTimeSelectionTable(functionObject, hydroGeoMetrics, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::hydroGeoMetrics::hydroGeoMetrics
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
fvMeshFunctionObject(name, runTime, dict),
// -- TODO add default values
debug_(dict.lookupOrDefault("debug",false)),
nBins_(dict.lookupOrDefault("nBins",label(20))),
nRadii_(dict.lookupOrDefault("nRadii",label(20))),
maxDist_(dict.lookupOrDefault("maxDist",scalar(1e10))),
fieldName_(dict.lookup("fieldName")),
operations_(dict.lookup("operations")),
vsfPtr_(),
vvfPtr_(),
vtfPtr_(),
vstfPtr_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::hydroGeoMetrics::~hydroGeoMetrics()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::hydroGeoMetrics::read(const dictionary& dict)
{
    return true;
}


bool Foam::functionObjects::hydroGeoMetrics::execute()
{
    return true;
}


bool Foam::functionObjects::hydroGeoMetrics::end()
{
    return true;
}


bool Foam::functionObjects::hydroGeoMetrics::write()
{
    bool sfM = computeMetrics(vsfPtr_);
    bool vfM = computeMetrics(vvfPtr_);
    bool tfM = computeMetrics(vtfPtr_);
    bool stfM = computeMetrics(vstfPtr_);

    if(sfM || vfM || tfM || stfM)
    {
        return true;
    }
    else
    {
        FatalErrorInFunction
        << "There is no field " << fieldName_ << endl
        << exit(FatalError);
    }

    return false;
}

//#include "twoPointCorrelations.C"
//#include "entropy.C"

// ************************************************************************* //
