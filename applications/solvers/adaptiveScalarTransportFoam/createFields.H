//################################################################//

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//################################################################//

Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading diffusivity DT\n" << endl;
dimensionedScalar DT
(
    transportProperties.lookup("DT")
);

Info<< "Reading outlet name\n" << endl;
const word outletName
(
    transportProperties.lookupOrDefault<word>
    ("outletBoundary", "outlet")
);

Info << "Reading porosity beta\n" << endl;
const dimensionedScalar beta
(
    transportProperties.lookupOrDefault<dimensionedScalar>
    ("beta", dimensionedScalar(dimless,scalar(1.)))
);

Info << "Reading controlDict\n" << endl;

scalar THRS = runTime.controlDict().lookupOrDefault<scalar>("THRS", 0.99);
scalar mTHRS = runTime.controlDict().lookupOrDefault<scalar>("mTHRS", 0.01);
word outletPatch = runTime.controlDict().lookupOrDefault<word>("outletPatch", "outlet");

//################################################################//

#include "createPhi.H"

#include "createFvOptions.H"

//################################################################//
