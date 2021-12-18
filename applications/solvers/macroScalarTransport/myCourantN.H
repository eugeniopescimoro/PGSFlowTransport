/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Global
    CourantNo

Description
    Calculates and outputs the mean and maximum Courant Numbers.

Developer
  Eugenio Pescimoro
\*---------------------------------------------------------------------------*/

scalar CoNum = 0.0;
scalar meanCoNum = 0.0;
scalar maxCoNumX = 0.0;
scalar xDim = 0.0;

// It computes the cell dimensions
const faceList & ff = mesh.faces();
const pointField & pp = mesh.points();

forAll ( mesh.C(), celli)
{
    const cell & cc = mesh.cells()[celli];
    labelList pLabels(cc.labels(ff));
    pointField pLocal(pLabels.size(), vector::zero);

    forAll (pLabels, pointi)
           pLocal[pointi] = pp[pLabels[pointi]];

    scalar xDim = Foam::max(pLocal & vector(1,0,0)) - Foam::min(pLocal & vector(1,0,0));
    scalar yDim = Foam::max(pLocal & vector(0,1,0)) - Foam::min(pLocal & vector(0,1,0));
    scalar zDim = Foam::max(pLocal & vector(0,0,1)) - Foam::min(pLocal & vector(0,0,1));
}

// maxCoNum is the maximum Courant number computed using velocity components and cell dimension. Others Courant computed as (flow rates)/(cell volumes) underestimate the maximum time step
{
    scalarField sumQ
    (
        fvc::surfaceSum(mag(q))().primitiveField()
    );

    CoNum = 0.5*gMax(sumQ/mesh.V().field())*runTime.deltaTValue();

    meanCoNum =
        0.5*(gSum(sumQ)/gSum(mesh.V().field()))*runTime.deltaTValue();

    // maxCoNumX = gMax(U.component(vector::X)/xDim)*runTime.deltaTValue();
}

Info<< "Courant Number mean: " << meanCoNum
    << " max: " << CoNum << endl;

// ************************************************************************* //
