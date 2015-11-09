// Nemo airfoil selection script 
//
// Diego Montero; Gustavo Violato; Fernando Valentini
// First release: Oct. 2015

clc
xdel(winsid())
clear

NEMO_DIR = getenv('NEMO_ROOT');
AIRFOIL_DB = NEMO_DIR + '/dbfiles/airfoils';

// Physical variables
// ==================
GRAV = 9.80665;   // gravity [m/s^2]

// Water temperature .vs. density LookUp-Table
// Temperature [degC]; Density [kg*m^-3]
// Source: Wikipedia
RHO_W_LUT = [10 999.7026;
	       15 999.1026;
	       20 998.2071;
	       22 997.7735;
	       25 997.0479;
	       30 995.6502;
	       40 992.2;
	       60 983.2]; 

// Water temperature .vs. dynamic viscosity LookUp-Table
// Temperature [degC]; Viscosity [mPa*s]
// Source: Wikipedia
MU_W_LUT = [10 1.308;
	      20 1.002;
	      30 0.7978;
	      40 0.6531;
	      50 0.5471;
	      60 0.4658];

// User fixed variables
// ====================

//Wing geometry
Sw  = 2.0^2/18; // Wing area  [m^2]
cr  = 159/1000; // Root chord [m]
ct  = 63/1000;  // Tip chord  [m]

w_temp = 21; // Water temperature [degC]
// Obs.  Temp varies annually between 17degC and 25degC @ Babitonga Bay
// Source: http://www.avesmarinhas.com.br/Mestrado%20Mario.pdf

v_cruise = 5.; // Cruise speed [m/s]

// Riders parameters
rider.mass = 85;          // [kg]
rider.pwr  = 250;         // [W]
rider.name = "Rider";

mass_hph = 15.; // Vehicle mass [kg]

// Auxiliary functions
// ===================

function rho_w = GetDens(temp)
  
  rho_w = interp1(RHO_W_LUT(:,1),RHO_W_LUT(:,2),temp);
  
endfunction

function visc_w = GetVisc(temp)
  
  visc_w = interp1(MU_W_LUT(:,1),MU_W_LUT(:,2),temp)*10^-3;
  
endfunction

// Calculations
// ============

// Rest of wing geometry
lbd = ct/cr;
mac = cr*2/3*(1+lbd+lbd^2)/(1+lbd);

Re_root = GetDens(w_temp)*v_cruise*cr/GetVisc(w_temp);
Re_mac  = Re_root*mac/cr;
Re_tip  = Re_root*ct/cr;

mass_tot = rider.mass + mass_hph;

CL_cruise = (mass_tot*GRAV)/(0.5*GetDens(w_temp)*v_cruise^2*Sw);

// Output
// ======
mprintf('Results:\n');
mprintf('----------------\n');
mprintf('Re@Root: %.0f\n',Re_root);
mprintf('Re@MAC:  %.0f\n',Re_mac);
mprintf('Re@Tip:  %.0f\n',Re_tip);
mprintf('Cruise CL: %.2f\n',CL_cruise);
mprintf('Re_MAC*sqrt(CL): %.2f\n',Re_mac*sqrt(CL_cruise));

// Change log
// ==========
// 2015-10-18 - Gustavo Violato
// ----------------------------
// First release
