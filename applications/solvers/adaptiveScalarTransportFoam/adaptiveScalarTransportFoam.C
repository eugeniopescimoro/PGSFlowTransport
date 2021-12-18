/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         |
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    adaptiveScalarTransportFoam

Description
    Solves the steady or transient transport equation for a passive scalar using
    an adaptive time step while setting a threshold as a lower limit for the
    passive scalar to be computed.

Developers
    Eugenio Pescimoro
    Matteo Icardi

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    //#include "createTimeControls.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"
    scalar maxCo = runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);

    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi))().primitiveField()
    );

    runTime.setDeltaT (maxCo/(0.5*gMax(sumPhi/mesh.V().field())));

    label patchId = mesh.boundaryMesh().findPatchID(outletPatch);
    const scalarField& TOut(T.boundaryField()[patchId]);
    const vectorField& Sf( mesh.Sf().boundaryField()[patchId]);
    const vectorField& Uout(U.boundaryField()[patchId]);

    scalarField fluxOut((Uout&Sf));
    //Info << "Flux out = " << gSum(fluxOut) << nl << endl;

    vector meanVel( (fvc::domainIntegrate(U)).value() / gSum(mesh.V()) );
    Info << "Mean vel = " << meanVel << nl << endl;

    scalar massFlag(0.0);

    while (simple.loop(runTime))
    {
        Info << "Adaptive time = " << runTime.timeName() << endl;

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(beta, T)
              + fvm::div(phi, T)
              - fvm::laplacian(beta*DT, T)
             ==
                fvOptions(T)
            );

            TEqn.relax();
            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T);
        }

        scalar totalMass( (fvc::domainIntegrate(T)).value() );
        scalar cFluxOut(gSum(TOut*(Sf&Uout))/gSum(fluxOut) );

        Info << "Total mass = " << totalMass << endl;
        Info << "Flux out = " << cFluxOut << nl << endl;

        if ( gSum(TOut*(Sf&Uout))/gSum(fluxOut) > THRS )
        {
            T.write();
            FatalErrorIn("adaptiveScalarTransportFoam")
                << "Concentration under the threshold" << endl
                << exit(FatalError);
        }

        scalar massFlagNew(fmod(totalMass,mTHRS));
        if ( massFlagNew < massFlag )
        {
          T.write();
        }
        massFlag = massFlagNew;

        runTime.write();
    }

    Info << "End\n" << nl << endl;

    return 0;
}


// ************************************************************************* //
