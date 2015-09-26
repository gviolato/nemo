// HPH preliminary design script
//
// Based on article "The 20-knot human-powered water craft" - Alec Brooks
// "Human Power" Magazine, Spring 1987, Vol6 N.1 - p.1
//
// Diego Montero; Gustavo Violato
// First release: Sep. 2015


// User defined variables
rho_w    = 1000;    // density of water 16C         [kg/m^3]
rho_a    = 1.225;  // density of air 15C           [kg/m^3]
CL_max_w = 1.1;    // max wing lift coeff          [-]
mass_tot = 85+15;  // total mass (rider + vehicle) [kg]
tc_w     = 0.13;   // wing thickness ratio         [-]
lambda_w = 0.4;    // wing taper ratio             [-]
eff_prop = 0.85;   // propeller efficiency         [-]
eff_mech = 0.95;   // mechanical drivetrain efficiency
f_i      = 1.48;   // induced drag factor          [-]
CD_0_w   = 0.008;  // wing profile drag coeff
CD_strut = 0.0085; // strut profile drag coeff
S_strut  = 0.0427; // strut area (submerged)       [m^2]
t_strut  = 0.028;  // strut thickness              [m]
CD_fw    = 0.009;  // front-wing drag coeff
S_fw     = 0.04;   // front-wing area              [m^2]
CD_spray = 0.24;   // spray drag coeff
CD_air   = 0.7;    // Drag coefficient of air-exposed comp.
S_air    = 0.6503; // Area of air-exposed comp.    [m^2]

// Design power level
Pwr_des  =  0.25;   // 30min from Whitt & Wilson    [kW]

// Fixed variables
GRAV = 9.80665;   // gravity [m/s^2]

// Drag Sum-Up (all sources except wing)

//     Interference                   Spray
Sref = (17*(tc_w)^2-0.05)*t_strut^2 + CD_spray*t_strut^2 + ...
       CD_air*S_air*rho_a/rho_w + CD_strut*S_strut + CD_fw*S_fw;

// As presented in the paper, we will build iso-lines of certain
// performance parameters (Max Speed, TO-Speed, Max wing-tip deflection)
// as a function of aspect-ratio and wing span

function A = MaxSpeed(speed_lvl, wing_span)
  
  e = 1/(1+0.015);
  q = 0.5*rho_w*(speed_lvl)^2;
  Di = f_i*(mass_tot*GRAV)^2/(%pi*e*q*wing_span^2);
  S = (Pwr_des*eff_prop*eff_mech/speed_lvl - q*Sref - Di)/(q*CD_0_w);
  A = wing_span^2/S;
  
endfunction

// Change log

// 2015-09-20 - Gustavo Violato
// First release

