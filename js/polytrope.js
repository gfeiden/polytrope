/* polytrope.js
   
   written by Gregory Feiden
*/

var Constants = {
    // Solar properties
    Msun: 1.989e33,      // solar mass in g
    Rsun: 6.956e10,      // solar radius in cm
    Lsun: 3.846e33,      // solar luminosity in erg/s
    
    // Physical constants
    m_proton: 1.672e-24, // mass of proton (g)
    G: 6.673e-8,         // Newton's gravitational G
    k_boltz: 1.381e-16,  // Boltzmann constant
    sigma: 5.671e-5      // Stefan-Boltzmann constant
}

function eos(Star, kind) {
    // define default EOS as ideal gas
    kind = kind || "ideal";

    if (kind == "degen") {
        var gamma = 4.0/3.0;
    }
    else {
        var gamma = 5.0/3.0;
    }
    
    // define for convenience
    var Mtot = Star.mass*Constants.Msun;
    var Rtot = Star.radius*Constants.Rsun;
    
    // central thermodynamic variables
    var M = [0.0];
    var rho = [-Mtot*Star.zero*Math.pow(Star.slope, -1)/(4.0*Math.PI*Math.pow(Rtot, 3))];
    var P = [Math.pow(Rtot/Star.zero, 2)*(4.0*Math.PI*Constants.G*Math.pow(rho[0], 2))/(Star.n + 1.0)];
    var K = P[0]/Math.pow(rho[0], gamma);
    var R = [Math.sqrt(((Star.n + 1.0)*P[0])/(4.0*Math.PI*Constants.G*Math.pow(rho[0], 2)))*Star.x[0]];

    if (kind === "ideal") {
        var T = [Constants.m_proton*P[0]/(Star.mu*Constants.k_boltz*rho[0])];
        for (i = 1; i < Star.y.length - 1; i++) {
            rho[i] = rho[0]*Math.pow(Star.y[i], Star.n);
            P[i] = K*Math.pow(rho[i], 1.0/Star.n + 1.0);
            T[i] = Constants.m_proton*P[i]/(Star.mu*Constants.k_boltz*rho[i]);
            R[i] = Math.sqrt(((Star.n + 1.0)*P[i])/(4.0*Math.PI*Constants.G*Math.pow(rho[i], 2)))*Star.x[i];
            M[i] = M[i - 1] + 4.0*Math.PI*Math.pow(R, 2)*rho[i]*(R[i] - R[i - 1]);
        }
    }
    else if (kind == "degen") {
        T = [0.0];
        for (i = 1; i < Star.y.length - 1; i++) {
            rho[i] = rho[0]*Math.pow(y[i], Star.n);
            P[i] = K*Math.pow(rho[i], 1.0/Star.n + 1.0);
            T[i] = 0.0;
            R[i] = Math.sqrt(((Star.n + 1.0)*P[i])/(4.0*Math.PI*Constants.G*Math.pow(rho[i], 2)))*Star.x[i];
            M[i] = M[i - 1] + 4.0*Math.PI*Math.pow(R, 2)*rho[i]*(R[i] - R[i - 1]);
        }
    }
    else {
        rho = ["Incorrect EOS"];
        P = ["Incorrect EOS"];
        T = ["Incorrect EOS"];
        R = ["Incorrect EOS"];
        M = ["Incorrect EOS"];
    }
    return [rho, P, T, R, M];
}

function zCalc(z, dz, h) {
    return z + h*dz/2.0;
}

function dzCalc(x, y, z, n) {
    return -Math.pow(y, n) - 2.0*z/x;
}

function solveRK4(x, y, z, dz, h, n) {
    // Solve the differential equations using RK order 4
    var i = 0;
    while (y[i] >= 0.0) {
        // DETERMINE K1, L1
        var k1 = h * z[i];
        var l1 = h * dzCalc(x[i], y[i], z[i], n);

        // DETERMINE K2, L2
        var k2 = h * zCalc(z[i] + l1/2.0, dz[i], h);
        var l2 = h * dzCalc(x[i] + h/2.0, y[i] + k1/2.0, z[i] + l1/2.0, n);

        // DETERMINE K3, L3
        var k3 = h * zCalc(z[i] + l2/2.0, dz[i], h);
        var l3 = h * dzCalc(x[i] + h/2.0, y[i] + k2/2.0, z[i] + l2/2.0, n);

        // DETERMINE K4, L4
        var k4 = h * zCalc(z[i] + l3, dz[i], h);
        var l4 = h * dzCalc(x[i] + h, y[i] + k3, z[i] + l3, n);

        // ITERATE ON X, Y, Z, I
        x[i + 1] = x[i] + h;
        y[i + 1] = y[i] + k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0;
        z[i + 1] = z[i] + l1/6.0 + l2/3.0 + l3/3.0 + l4/6.0;
        dz[i + 1] = dzCalc(x[i], y[i], z[i], n);

        i++;
    }

    return [x, y, z, dz];
}

function starInit(mass, teff, luminosity, n, X, Y, Z) {
    // Initialize star object
    var Star = {
        mass: mass,
        teff: teff,
        luminosity: luminosity,
        radius: Math.sqrt(luminosity*Constants.Lsun/(4.0*Math.PI*Constants.sigma*Math.pow(teff, 4)))/Constants.Rsun,
        n: n,
        X: X,
        Y: Y,
        Z: Z,
        mu: 1.0/(2.0*X + 0.75*Y + 0.5*Z)
    }
    Star.rk4_x = [0.001];
    Star.rk4_y = [1.0 - Math.pow(Star.rk4_x[0], 2)/6.0];
    Star.rk4_z = [-Star.rk4_x[0]/3.0];
    Star.rk4_dz = [dzCalc(Star.rk4_x[0], Star.rk4_y[0], Star.rk4_z[0], Star.n)];

    return Star;
}

function polytrope(mass, teff, luminosity, n, Y, Z, h, eosKind) {
    // check that the proper number of arguments were sent
    
    // define default values
    mass = mass || 1.0;
    teff = teff || 5776.0;
    luminosity = luminosity || 1.0;
    n = n || 1.0;
    Y = Y || 0.28;
    Z = Z || 0.02;
    h = h || 0.001;
    eosKind = eosKind || "ideal";
    
    var Star = starInit(mass, teff, luminosity, n, Y, Z);
    var solution = solveRK4(Star.rk4_x, Star.rk4_y, Star.rk4_z, Star.rk4_dz, h, n);
    Star.x = solution[0];
    Star.y = solution[1];
    Star.z = solution[2];
    Star.dz = solution[3];
    
    // linear extrapolate to find zero point
    var len = Star.y.length - 1;
    Star.slope = (Star.y[len - 2] - Star.y[len - 1])/(Star.x[len - 2] - Star.x[len - 1]);
    Star.zero = Star.x[len - 1] - Star.y[len - 1]/Star.slope;

    var props = eos(Star, eosKind);
    Star.rho = props[0];
    Star.P = props[1];
    Star.T = props[2];
    Star.R = props[3];
    Star.M = props[4];
    
    return Star;
}

