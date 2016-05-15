// HPH preliminary design script
//
// Based on article "The 20-knot human-powered water craft" - Alec Brooks
// "Human Power" Magazine, Spring 1987, Vol6 N.1 - p.1
//
// Diego Montero; Gustavo Violato; Fernando Valentini
// First release: Sep. 2015

clc
xdel(winsid())
clear

// Script flags
DO_CALC=1;
DO_PLOTS=1;

// Fixed variables
GRAV = 9.80665;   // gravity [m/s^2]
FT2M = 0.3048;    // feet to meter

// Riders parameters
gustavo.mass = 87;
gustavo.pwr  = 215;
gustavo.name = "Gustavo";

fernando.mass = 88;
fernando.pwr  = 250;
fernando.name = "Valentini";

diego.mass = 72;
diego.pwr = 275;
diego.name = "Diego";

// Chosen wing dimensions
prj_span = 2.0;    // [m]
prj_aspect = 18;    // [m]

// Rider
rider = gustavo;

// User defined variables
rho_w    = 1000;    // density of water 16C        [kg/m^3]
rho_a    = 1.225;  // density of air 15C           [kg/m^3]
CL_max_w = 1.2;    // max wing lift coeff          [-]
mass_hph = 15;     // vehicle mass                 [kg]
tc_w     = 0.13;   // wing thickness ratio         [-]
lambda_w = 0.4;    // wing taper ratio             [-]
eff_prop = 0.85;   // propeller efficiency         [-]
eff_mech = 0.95;   // mechanical drivetrain efficiency
f_i      = 1.48;   // induced drag factor          [-]
CD_0_w   = 0.007;  // wing profile drag coeff
CD_strut = 0.0085; // strut profile drag coeff
S_strut  = 0.0427; // strut area (submerged)       [m^2]
t_strut  = 0.028;  // strut thickness              [m]
CD_fw    = 0.009;  // front-wing drag coeff
S_fw     = 0.04;   // front-wing area              [m^2]
CD_spray = 0.24;   // spray drag coeff
CD_air   = 1.0;    // Drag coefficient of air-exposed comp.
S_air    = 0.35; // Area of air-exposed comp.    [m^2]
// "20-Knot" article values are CD_air=0.7 and S_air=0.6503
// Recumbent "Robinho" area is 0.377m^2
// It seems that a 'drag area" of 0.35m^2 for recumbent is reasonable

// Drag Sum-Up (all sources except wing)

Drag_q.Interference = (17*(tc_w)^2-0.05)*t_strut^2;
Drag_q.Spray = CD_spray*t_strut^2;
Drag_q.Air   = CD_air*S_air*rho_a/rho_w;
Drag_q.Strut = CD_strut*S_strut;
Drag_q.Front = CD_fw*S_fw;

Sref = 0;
fn = fieldnames(Drag_q);
for n=1:length(length(fn))
  Sref = Sref + Drag_q(fn(n));
end

// Deflection data (Taken from Figure 6.5, chapter 6 of "Human-Powered
// Vehicles", p. 87)
// 0.08 tip deflection for Fiberglass
DEFL = [2 (1.5/7.3)*30+5;
	(3.75/11.4)*7+2 (2.7/7.3)*30+5;
	9 (4.4/7.3)*30+5];

// As presented in the paper, we will build iso-lines of certain
// performance parameters (Max Speed, TO-Speed, Max wing-tip deflection)
// as a function of aspect-ratio and wing span

function Pwr = OpPointPwr(aspect,wing_span,speed,rider)
  
  mass_tot = rider.mass + mass_hph;
  
  e = 1/(1+0.015);
  q = 0.5*rho_w*speed^2;
  S = wing_span^2/aspect;
  Di = f_i*(mass_tot*GRAV)^2/(%pi*e*q*wing_span^2);
  S_drag = CD_0_w*S+Sref;
  
  Pwr = speed*(q*S_drag + Di)/(eff_prop*eff_mech);
  
endfunction

function V = MaxSpeed(aspect,wing_span,rider)
  
  mass_tot = rider.mass + mass_hph;
  Pwr_des  = rider.pwr;
  
function y=dragequation(x)
  
  v = x(1);
  
  e = 1/(1+0.015);
  q = 0.5*rho_w*v^2;
  S = wing_span^2/aspect;
  Di = f_i*(mass_tot*GRAV)^2/(%pi*e*q*wing_span^2);
  
  y = Pwr_des*eff_prop*eff_mech - v*q*(CD_0_w*S+Sref) - v*Di;
  
endfunction

  [xres,val,info] = fsolve([8],dragequation);
  if info==1
    V = xres(1);
  else
    V = %nan;
  end
  
endfunction

function [Pwr, Vel] = MinPower(aspect,wing_span,rider)

  mass_tot = rider.mass + mass_hph;
  
function y = minpowersystem(x)
  
  y = zeros(2,1);
  
  v   = x(1);
  pwr = x(2);

  e = 1/(1+0.015);
  q = 0.5*rho_w*v^2;
  S = wing_span^2/aspect;
  Di = f_i*(mass_tot*GRAV)^2/(%pi*e*q*wing_span^2);
  S_drag = CD_0_w*S+Sref;
  
  y(1) = pwr*eff_prop*eff_mech - v*(q*S_drag+Di);
  y(2) = 3*q*(S_drag)- Di;
  
endfunction

  [xres,val,info] = fsolve([3;1.5],minpowersystem);
  if info==1
    Vel = xres(1);
    Pwr = xres(2);
  else
    Vel = %nan;
    Pwr = %nan;
  end

endfunction //MinPower

function V = TOSpeed(aspect,wing_span,rider)
  
  mass_tot = rider.mass + mass_hph;
  V = sqrt(aspect*mass_tot*GRAV/...
	   (0.5*rho_w*0.9*CL_max_w*wing_span^2));  
endfunction // TOSpeed

//Axiliary functions

function map = mapspace(func,xvec,yvec,extras)
  
  rider = extras(1).entries;
  
  xi=1;
  for x=xvec
    yi=1;
    for y=yvec
      map(yi,xi)=func(x,y,rider);
      yi=yi+1;
    end
    xi=xi+1;
  end
  
endfunction

// Plotting and printing functions
function setmyaxprops()
  
  ax = gca();
  
  myfontstyle = 5;
  
  ax.filled = "on";
  ax.background = color('white');

  ax.grid = [color('gray') ,color('gray')];
  ax.grid_style = [7,7]; // short-dashed lines
  
  ax.font_style = myfontstyle;
  ax.font_size = 2;
  
  ax.y_label.font_style = myfontstyle;
  ax.y_label.font_size = 3;
  ax.x_label.font_style = myfontstyle;
  ax.x_label.font_size = 3;
  
  ax.title.font_style = myfontstyle;
  ax.title.font_size = 4;
  
endfunction // setmyaxprops

function cr = PlotWing(aspect,wing_span,lambda)
  
  S  = wing_span^2/aspect;
  cr = 2*S/(wing_span*(1+lambda));
  
  x_tip_l  = -wing_span/2;
  x_tip_r  = wing_span/2;
  
  y_tip_lw = -lambda*cr/2;
  y_tip_up = lambda*cr/2;
  
  pol_x = [x_tip_l x_tip_l 0 x_tip_r x_tip_r 0];
  pol_y = [y_tip_lw y_tip_up cr/2 y_tip_up y_tip_lw -cr/2];
  
  figure();
  xpoly(pol_x,pol_y,"lines",1);
  xfpolys(pol_x',pol_y',color('white'));
  ax = gca();
  ax.isoview = "on";
  
endfunction // PlotWing

function PlotDragPie(drags)
  
  x   = [];
  sp  = [];
  txt = [];
  fn = fieldnames(drags);
  for n=1:length(length(fn))
    x = [x drags(fn(n))];
    txt = [txt fn(n)];
    if length(strstr(fn(n),'Wing'))>0
      sp = [sp 1];
    else
      sp = [sp 0];
    end
  end
  figure()
  pie(x,sp,txt);
  
endfunction

// Script calculations

if DO_CALC

  // Extra variables for mapping functions
  extras = cell(1,1);
  extras(1).entries = rider;
  
  // First define reasonable intervals for aspect ratio and
  // wing span
  aspects = linspace(5,35,75);
  span_interval = linspace(1.05,2.3,75);
  
  // Calculate isolines of the interesting performance parameters
  // as function of wing dimensions using the above-defined functions
  V_max   = mapspace(MaxSpeed,aspects,span_interval,extras);
  Pwr_min = mapspace(MinPower,aspects,span_interval,extras);
  V_TO    = mapspace(TOSpeed,aspects,span_interval,extras);
  
end

// Plot results
if DO_PLOTS
  
  speeds    = [4.5 5 5.5 6];
  minpwrlvl = [125 150 175 200];
  speedsTO  = [2.5 2.75 3 3.25];

  N_speeds = length(speeds);
  N_pwr    = length(minpwrlvl);
  N_TO     = length(speedsTO);
  
  // Isoline plots
  
  figure();
  
  contour2d(span_interval,aspects,V_max,speeds,...
	    style=color('blue')*ones(1,N_speeds));
  contour2d(span_interval,aspects,Pwr_min,minpwrlvl,...
	    style=color('red')*ones(1,N_pwr));
  contour2d(span_interval,aspects,V_TO,speedsTO,...
	    style=color('green')*ones(1,N_TO));
  plot2d(DEFL(:,1)*FT2M,DEFL(:,2),color('magenta'));
  plot2d(prj_span,prj_aspect,-3);
  
  xgrid();
  xlabel('Wing Span [m]');
  ylabel('Aspect Ratio');
  xtitle('Performance Curves. Rider: '+rider.name);
  
  legends(['Take-off speed','Cruise speed','Min power','Tip deflection'], [color('green') color('blue') color('red') color('magenta')], 'ur');

  setmyaxprops();
  ax = gca();
  ax.tight_limits = 'on';
  ax.data_bounds = [1.1 5; 2.3 30];

  // Wing Plot
  cr = PlotWing(prj_aspect,prj_span,lambda_w);
  
  // Drag buil-up @ 5m/s
  mass_tot = rider.mass + mass_hph;
  e = 1/(1+0.025);
  q = 0.5*rho_w*5^2;
  Drag_q.WingCD0 = CD_0_w*prj_span^2/prj_aspect;
  Drag_q.WingDi  = f_i*(mass_tot*GRAV)^2/(%pi*e*q^2*prj_span^2);
  PlotDragPie(Drag_q);
  
end

// Compute results
to_speed = TOSpeed(prj_aspect,prj_span,rider);
to_power = OpPointPwr(prj_aspect,prj_span,to_speed,rider);
cruise_speed = MaxSpeed(prj_aspect,prj_span,rider);
old_pwr = rider.pwr;
rider.pwr = 750;
max_speed = MaxSpeed(prj_aspect,prj_span,rider);
rider.pwr = old_pwr;
[min_power minpwrvel] = MinPower(prj_aspect,prj_span,rider);

// Print results
mprintf('Results:\n');
mprintf('----------------\n');
mprintf('Chosen wind span[m]/aspect ratio[-]: %.2f, %.1f\n\n',prj_span,prj_aspect);
mprintf('-Rider----------\n');
mprintf('Name     : %s\n',rider.name);
mprintf('Mass [kg]: %.1f\n',rider.mass);
mprintf('Power [W]: %.1f\n\n',rider.pwr);
mprintf('-Wing-----------\n');
mprintf('Wing area [m^2]: %.3f\n',prj_span^2/prj_aspect);
mprintf('Root Chord  [m]: %.3f\n',cr);
mprintf('Tip Chord   [m]: %.3f\n\n',lambda_w*cr);
mprintf('-Performance----\n');
mprintf('Take-off speed [m/s]: %.1f\n',to_speed);
mprintf('Take-off power   [W]: %.1f\n',to_power);
mprintf('Cruise Speed   [m/s]: %.1f\n',cruise_speed);
mprintf('Range           [km]: %.2f\n',cruise_speed*1.0*3600/1000);
mprintf('Speed@750W     [m/s]: %.1f\n',max_speed);
mprintf('Min power        [W]: %.1f\n',min_power);
mprintf('Speed@MinPwr   [m/s]: %.1f\n',minpwrvel);


// Change log

// 2015-10-17 - Diego Montero
// ----------------------------
// Added legends to the line functions :)

// 2015-09-29 - Gustavo Violato
// ----------------------------
// Change original line functions for height-maps of desired variables
// Plot of wing planform and drag build-up for design point as a pie chart

// 2015-09-20 - Gustavo Violato
// ----------------------------
// First release

