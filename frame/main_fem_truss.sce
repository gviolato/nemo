// Nemo Frame preliminary design script
//
// Finite Element Aproach for frame design
//.Working Beam elements for 2D and 3D with Mesh refinement
//
// Diego Montero; Fernando Valentini; Gustavo Violato; 
// First release: Oct. 2015 
clc; xdel(winsid()); clear;

//Add the NEMO_ROOT enviroment variable as the project folder in your computer to run the codes
nemo_root = getenv('NEMO_ROOT');
cd(nemo_root+'\frame');
exec(pwd()+'\macros\library.sce');

//Problem choice + Properties
exec('2D_10-bar-truss.sce');

//Mesh Refinement, input 0 to skip
Ref=6000;
[Ref_conec,Ref_coord,Ref_Prop]=refine(conec,coord,Prop,Ref);

//Solve truss
[Results,Udef]=solve_truss(Ref_conec,Ref_coord,Ref_Prop,F);

//Display Results
printf("  Node   Stress(MPa) Lin.Disp.(mm)")
disp(Results);

//Plot results
plottruss(Ref_conec,Ref_coord,'green');
Def_coord=Ref_coord+10*Udef;
plottruss(Ref_conec,Def_coord,'red');
legends(['before loading','after loading'], [color('green') color('red')], 'ur');

// Change log

//2015-10-17 - Diego Montero
// ----------------------------
//Fixed a potential bug on the DOF restrictions;
//Code already benchmarked against other FEA softwares, everything is working fine so far;

//2015-10-14 - Diego Montero
// ----------------------------
//Implemented Mesh refinement;

//2015-10-10 - Diego Montero
// ----------------------------
//Implemented 3D function;
//Program splited into main, problem, plot and solve files

// 2015-10-04 - Diego Montero
// ----------------------------
// First release, working for 2D truss with simple plot in a single file.
