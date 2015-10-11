//2D Triangle Problem

//Material properties & Geometry
Em=200e3; //Young Modulus in MPa (1020 Steel = 200e3)
Im=489.86; //Moment of Inertia in mm4 (10mm diameter bar = 489.86)

 //Area in mm2 (10mm diameter bar = 78.4787)
Am=78.4787;
//If any properties change for each element, use matrix input
Gm=0;
Jm=0;
//Nodes coordinates, each line is a node n located at x & y
coord = [   0 0
20 0
40 0
60 0
80 0
100 0
];

//Conectivity Matrix, each line is an element connecting nodes a & b
conec = [1 2
2 3
3 4
4 5
5 6
];

//Forces, applied to matching DOF
    F(17)=-100;

//Restricted DOFs (also for boundary conditions)
glrest=[1; 2; 3];
