//2D_10-bar-truss

//Nodes coordinates, each line is a node n located at x & y
coord = [18288 9144 
18288 0
9144 9144
9144 0
0 9144
0 0 ];

//Conectivity Matrix, each line is an element connecting nodes a & b
conec = [5 3
3 1
6 4
4 2
4 3
2 1
5 4
6 3
3 2
4 1
];

//Forces, applied to matching DOF
    F(5)=-444822.16;
    F(11)=-444822.16;

//Restricted DOFs (also for boundary conditions)
glrest=[13,14,16,17];//7*3,8*3,9*3,10*3,11*3,12*3,13*3,14*3,15*3,16*3];
Am=[1000];
d=sqrt(4*Am/%pi) ;
Im=(%pi*d^4) / 64;
Em=68.95E3;
Gm=0;
Jm=0;


//If any properties change for each element, use matrix input

//Create Properties Matrix
[Prop]=Properties(Am,Em,Im,Gm,Jm,conec);
