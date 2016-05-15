// HPH preliminary hull shape
//
// Some parametric hull shape tests
//
// Diego Montero; Gustavo Violato; Fernando Valentini
// First release: Nov. 2015

clc
xdel(winsid())
clear

// Script flags
PLOT = 1;

// User defined variables

L   = 2.00;    // hull length in m
Prof = 0.20;
Aber = 0.25/2;

C_L = Prof/L;
B_L = Aber/L;

NS  = 15;

// Auxiliary functions

function p = normalized_pol2(x,frac)
  
  p = 4*frac*(x.^2 - 1/4);
  
endfunction

function p = normalized_pol4(x,frac)
// normalized_pol: 
  
// x -> number between -1/2 and 1/2: normalized length
// frac -> fraction of length correspondent to maximum width
  
  p = x.^4 + (4*frac-1/4)*x.^2 - frac;
  
endfunction

function l = line_coord(s,len,nf)
  
  l = len*nf( (s-len/2)/len );
  
endfunction

function c = cross_section(s,len,nf_c,nf_b,cross_fun)
  
  xN = (s-len/2)/len;
  
  b = len*nf_b(xN);
  c = len*nf_c(xN);
  
  bN = 0:0.01:1;
  y = c*cross_fun(bN);
  
  c = [bN*b; y];
  
endfunction

function d = displacement(len,nf_c,nf_b,area_fun)

  s = linspace(0,len);
  b = line_coord(s,len,nf_b);
  c = line_coord(s,len,nf_c);
  a = area_fun(b,c);
  
  d = intsplin(s,a);
  
endfunction
  
// Script calculations

// Points defining the stations
s = 0:0.01:L;

// Keel coordinates
deff('[pk]=keel_fun(x)','pk=normalized_pol4(x,C_L)');
keel_y = 0*s;
keel_z = line_coord(s,L,keel_fun);

// Water line coordinates
deff('[pw]=wl_fun(x)','pw=normalized_pol2(x,B_L)');
wl1_y = line_coord(s,L,wl_fun);
wl1_z = 0*s;

wl2_y = -1*wl1_y;
wl2_z = 0*s;


// Cross sections
deff('[c]=cross_fun(x)','c=sqrt(1-x.^2)');
deff('[a]=area_fun(b,c)','a=%pi*b.*c/2');
sts = linspace(0,L,NS+2);
sts = sts(2:$-1);

css = hypermat([2,201,NS]);
i = 1;
for cs=sts
  c  = cross_section(cs,L,keel_fun,wl_fun,cross_fun);
  ycs = [-1*c(1,$:-1:2) c(1,:)];
  zcs = [c(2,$:-1:2) c(2,:)];
  css(:,:,i) = [ycs; zcs];
  i = i + 1;
end

if PLOT
  // Plotting
  param3d(s,keel_y,keel_z);
  param3d(s,wl1_y,wl1_z);
  param3d(s,wl2_y,wl2_z);
  for si=1:NS
    x_cs = sts(si)*ones(1,201);
    param3d(x_cs,css(1,:,si),css(2,:,si));
  end
  a = gca();
  a.isoview = "on";
end


mprintf('Results:\n');
mprintf('---------------\n');
mprintf('Displacement: %.3f\n',displacement(L,keel_fun,wl_fun,area_fun));
mprintf('Width (Boca): %.3f\n',2*B_L*L);
mprintf('Depth (Pontal): %.3f\n',C_L*L);
