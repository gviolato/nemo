//3D Pyramid Problem

//Nodes coordinates, each line is a node n located at x, y & z
coord = [   0 0 0
 18288 0 0
18288 18288 0 
 0 18288 0
9144 9144 18288
];

//Conectivity Matrix, each line is an element connecting nodes a & b
conec = [1 2
        2 3
        3 4
        4 1
        5 1
        5 2
        5 3
        5 4
];

//Forces, applied to matching DOF
    F(27)=-444822.16 ; //in this case F in y of node 5 (DOF 25)

//Restricted DOFs
glrest=[1; 2; 3;4;5;6;9; 15;21];

//Material properties & Geometry
po=0.35;
Am=[1000];
d=sqrt(4*Am/%pi) ;
Im=(%pi*d^4) / 64;
Em=68.95E3;
Gm=Em/(2*(1+po));
Jm=2*Im;

//If any properties change for each element, use matrix input

//Create Properties Matrix
[Prop]=Properties(Am,Em,Im,Gm,Jm,conec);


