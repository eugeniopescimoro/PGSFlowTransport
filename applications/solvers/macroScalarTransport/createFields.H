//################################################################//

Info<< "Reading field c\n" << endl;
volScalarField c
(
    IOobject
    (
        "c",
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

Info<< "Reading diffusivity Dm\n" << endl;
dimensionedScalar Dm
(
    transportProperties.lookup("Dm")
);

Info<< "Reading outlet name\n" << endl;
const word outletName
(
    transportProperties.lookupOrDefault<word>
    ("outletBoundary", "outlet")
);

Info << "Reading porosity phi\n" << endl;
const dimensionedScalar phi
(
    transportProperties.lookupOrDefault<dimensionedScalar>
    ("phi", dimensionedScalar(dimless,scalar(1.)))
);

Info << "Reading mechanical disp coeff alpha\n" << endl;
const dimensionedScalar alpha
(
    transportProperties.lookupOrDefault<dimensionedScalar>
    ("alpha", dimensionedScalar(dimless,scalar(1.)))
);

Info << "Reading controlDict\n" << endl;
scalar THRS = transportProperties.lookupOrDefault<scalar>("THRS", 0.99);
scalar mTHRS = transportProperties.lookupOrDefault<scalar>("mTHRS", 0.01);
word outletPatch = runTime.controlDict().lookupOrDefault<word>("outletPatch", "outlet");

Info<< "Reading/calculating face flux field q\n" << endl;

surfaceScalarField q
(
    IOobject
    (
        "q",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U)
);

//################################################################//

//#include "createPhi.H"

#include "createFvOptions.H"

//################################################################//
