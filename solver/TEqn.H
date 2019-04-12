{
    alphat = turbulence->nut()/Scht;
    alphat.correctBoundaryConditions();

    volScalarField alphaEff("alphaEff", turbulence->nu()/Sch + alphat);

    // update update the particle velocity
    U_p = U-(1.0-beta)*T_p*(fvc::div(fvc::flux(U), U)-g); // Dr. bala's formula
    U_p.correctBoundaryConditions();

    Info<< "Calculating particle face flux field phi_p\n" << endl;
    surfaceScalarField phi_p
    (
        "phi_p",
        fvc::flux(U_p)
    );

    fvScalarMatrix TEqn
    (
        fvm::div(phi_p, T)
      - fvm::laplacian(alphaEff, T)
     ==
        fvOptions(T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);
    rhok = 1.0 + T*(rho_p/rho_f -1);
}
