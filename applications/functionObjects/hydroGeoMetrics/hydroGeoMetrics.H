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

Class
    Foam::functionObjects::fieldMetrics

Description
    This function object compute spatial statistics of a random field

    Example of function object specification:
    \verbatim
    fieldMetrics1
    {
        type           fieldMetrics;
        libs ("libfieldMetricsFunctionObject.so");

        nr             100;
        Lr             1.0;
        fieldName      K;
    }
    \endverbatim

Usage
    \table
        Property     | Description                | Required | Default value
        type         | type name: fieldMetrics    | yes      |
        nr           | number of radial samples   | yes      |
        Lr           | Max radial distnace        | yes      |
        fieldName    | Name of random field       | yes      |
    \endtable

SourceFiles
    fieldMetrics.C

\*---------------------------------------------------------------------------*/

#ifndef hydroGeoMetrics_H
#define hydroGeoMetrics_H

#include "fvMeshFunctionObject.H"
#include "fvCFD.H"
#include "graph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class fieldMetrics Declaration
\*---------------------------------------------------------------------------*/

class hydroGeoMetrics
:
    public fvMeshFunctionObject
{
    // Private Data

    //- Debug
    Switch debug_;

    //- Number of bins
    label nBins_;

    //- Number of radial bins
    label nRadii_;

    //- Maximum radial distance
    scalar maxDist_;

    //- Name of random field
    word fieldName_;

    //- Operaations to perform
    wordList operations_;

    //- Store fields if required
    autoPtr<volScalarField> vsfPtr_;
    autoPtr<volVectorField> vvfPtr_;
    autoPtr<volTensorField> vtfPtr_;
    autoPtr<volSymmTensorField> vstfPtr_;

    template<class Type>
    bool checkField
    (
        autoPtr<GeometricField<Type, fvPatchField, volMesh>>& fPtr
    );

    template<class Type>
    bool computeMetrics
    (
        autoPtr<GeometricField<Type, fvPatchField, volMesh>>& fPtr
    );

    template<class Type>
    void meanAndVariance(const GeometricField<Type, fvPatchField, volMesh>& K);

    void twoPointCorrelation(const volScalarField& K);
    void twoPointCorrelation(const volVectorField& K);
    void twoPointCorrelation(const volSymmTensorField& K);
    void twoPointCorrelation(const volTensorField& K);

    void entropy(const volScalarField& K);
    void entropy(const volVectorField& K);
    void entropy(const volSymmTensorField& K);
    void entropy(const volTensorField& K);

    template<class Type>
    void meanVel(const GeometricField<Type, fvPatchField, volMesh>& U);

    void velocity(const volVectorField& U);

    

public:

    //- Runtime type information
    TypeName("hydroGeoMetrics");


    // Constructors

        //- Construct from Time and dictionary
        hydroGeoMetrics
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        hydroGeoMetrics(const hydroGeoMetrics&) = delete;


    //- Destructor
    virtual ~hydroGeoMetrics();


    // Member Functions

        //- Read the fieldMetrics data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Write the fieldMetrics
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const hydroGeoMetrics&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
  #include "hydroGeoMetricsTemplates.C"
  #include "meanAndVariance.C"
  #include "meanVel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
