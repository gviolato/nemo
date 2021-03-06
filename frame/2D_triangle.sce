//2D Triangle Problem

//Nodes coordinates, each line is a node n located at x & y
coord = [   0 0
            100 0
            50 75
];

//Conectivity Matrix, each line is an element connecting nodes a & b
conec = [1 2
2 3
3 1
];

//Forces, applied to matching DOF
F(8)=-30000; //in this case F in y of node 3 (DOF 8)

//Restricted DOFs (also for boundary conditions)
glrest=[1; 2; 5];

//Material properties & Geometry
Em=200e3; //Young Modulus in MPa (1020 Steel = 200e3)
Im=489.86; //Moment of Inertia in mm4 (10mm diameter bar = 489.86)
Gm=Em;
Jm=2*Im;
Am=[78.4787]; //Area in mm2 (10mm diameter bar = 78.4787)

//If any properties change for each element, use matrix input

//Create Properties Matrix
[Prop]=Properties(Am,Em,Im,Gm,Jm,conec);
