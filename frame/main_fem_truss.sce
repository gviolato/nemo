// Nemo Frame preliminary design script
//
// Finite Element Aproach for frame design
//.Working Beam elements for 2D and 3D with Mesh refinement
//
// Diego Montero; Fernando Valentini; Gustavo Violato; 
// First release: Oct. 2015 
clc; xdel(winsid()); clear;

//Change directory and call functions
cd('C:\Users\dimon\Dropbox\GitHub\nemo\frame');
exec(pwd()+'\macros\library.sce');

//Problem choice + Properties
exec('3D_tetrahedron.sce');

//Mesh Refinement, input 0 to skip
Ref=10;
[Ref_conec,Ref_coord,Ref_Prop]=refine(conec,coord,Prop,Ref);

//Solve truss
[Results,Udef]=solve_truss(Ref_conec,Ref_coord,Ref_Prop,F);

//Display Results
printf("  Node   Stress(MPa) Lin.Disp.(mm)")
disp(Results);

//Plot results
plottruss(Ref_conec,Ref_coord,'green');
Def_coord=Ref_coord+10*Ref*Udef;
plottruss(Ref_conec,Def_coord,'red');
legends(['before loading','after loading'], [color('green') color('red')], 'ur');

// Change log

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
