singlePhaseTransportModel laminarTransport(U, phi);

// particle density
dimensionedScalar rho_p("rho_p", dimDensity, laminarTransport);
// fluid density
dimensionedScalar rho_f("rho_f", dimDensity, laminarTransport);
// dummy fluid viscocity 
dimensionedScalar nu_p("nu_p", laminarTransport.lookup("nu"));
// Laminar Schmidt number
dimensionedScalar Sch("Sch", dimless, laminarTransport);
// Turbulent Schmidt number
dimensionedScalar Scht("Scht", dimless, laminarTransport);
// particle diameter
dimensionedScalar d_p("d_p", dimLength, laminarTransport);
