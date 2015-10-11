// Nemo Frame preliminary design script
//
// Finite Element Aproach for frame design
//.Working Beam elements for 2D and 3D
//
// Diego Montero; Fernando Valentini; Gustavo Violato; 
// First release: Oct. 2015 
clc; xdel(winsid()); clear;

//Change directory and call funcs
cd('C:\Users\dimon\Dropbox\GitHub\nemo\frame');
exec(pwd()+'\macros\library.sce');

//Problem choice
exec('2D_beam.sce');

[Results,Udef]=solve_truss(conec,coord,Em,Im,Am,Gm,Jm,F);

//Display Results
printf("  Node   Stress(MPa) Lin.Disp.(mm)")
disp(Results);

//Plot results
plottruss(conec,coord,'green');
coord=coord+10*Udef;
plottruss(conec,coord,'red');
legends(['before loading','after loading'], [color('green') color('red')], 'ur');

// Change log

//2015-10-10 - Diego Montero
// ----------------------------
//Implemented 3D function;
//Program splited into main, problem, plot and solve files

// 2015-10-04 - Diego Montero
// ----------------------------
// First release, working for 2D truss with simple plot in a single file.
