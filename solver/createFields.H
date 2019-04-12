Info<< "Reading physical properties\n" << endl;

Info<< "Reading particle volumetric concentration field T\n" << endl;
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

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
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

#include "createPhi.H"

#include "readTransportProperties.H"

Info<< "Creating turbulence model\n" << endl;
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

// Kinematic density for buoyancy force/ or the density ratio rho/rho_0
volScalarField rhok
(
    IOobject
    (
        "rhok",
        runTime.timeName(),
        mesh
    ),
    1.0 + T*(rho_p/rho_f -1)
);


// kinematic turbulent diffusivity m2/s
Info<< "Reading field alphat\n" << endl;
volScalarField alphat
(
    IOobject
    (
        "alphat",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

///////////////////////////////////////
dimensionedScalar& nu_ = nu_p;

//Calculate the particle terminal velocity
scalar U_t_mag = mag(g.value())*Foam::pow(d_p.value(),2)*(rho_p.value()-rho_f.value())/(18*nu_.value()*rho_f.value()); // Stokes Law
scalar Re_p = U_t_mag*d_p.value()/nu_.value(); // Calculate the particle Reynolds number
scalar C_d; // drag coefficient
scalar T_p_mag = Foam::pow(d_p.value(),2)*rho_p.value()/(18*nu_.value()*rho_f.value()); // Stokes particle time scale magnitude

Info<< "U_t_mag = " << U_t_mag << " T_p_mag = "  << T_p_mag << endl;

if (Re_p > 1.0)
{
	Info<< "Schiller and Naumann drag correction!!!\n" << endl;
 	for (int i=0; i<50; i++)
 	{
		C_d = 24.0*(1.0+0.15*Foam::pow(Re_p,0.687))/Re_p; // Schiller and Naumann finite reynolds drag correction
		C_d = max(C_d, 0.44); // limites the drag coefficent, this should valid until the drag crisis.
		U_t_mag = Foam::pow((4.0*mag(g.value())*d_p.value()/3.0/C_d*((rho_p.value()-rho_f.value())/rho_f.value())),0.5); //update the terminval velocity
		Re_p = U_t_mag*d_p.value()/nu_.value(); // Calculate the particle Reynolds number
		T_p_mag = Foam::pow(d_p.value(), 2.0)*(rho_p.value()/rho_f.value())/(18.0*nu_.value()*(1.0+0.15*Foam::pow(Re_p,0.687)));
		Info<< "U_t_mag = " << U_t_mag << " T_p_mag = "  << T_p_mag << endl;
 	}
}

dimensionedScalar T_p 
(
"T_p",
dimensionSet(0, 0, 1, 0, 0, 0, 0),
T_p_mag
);

scalar beta = 3.0/(2.0*(rho_p.value()/rho_f.value())+1.0);

Info<< "Create particle velocity feild U_p" << " \n" << endl;
volVectorField U_p
(
    IOobject
    (
        "U_p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
//////////////////////////////////////////////

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rhok*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    simple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

mesh.setFluxRequired(p_rgh.name());

#include "createMRF.H"
