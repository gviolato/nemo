//3D Triangle Problem

//Nodes coordinates, each line is a node n located at x, y & z
coord = [   0 0 0
            -70.71068 0 -70.71068
            -35.35534  75 -35.35534
];

//Conectivity Matrix, each line is an element connecting nodes a & b
conec = [1 2
2 3
1 3
];

//Forces, applied to matching DOF
F(14)=-30000;

//Restricted DOFs (also for boundary conditions)
glrest=[1; 2; 3; 8];

//Material properties & Geometry
Em=200e3; //Young Modulus in MPa (1020 Steel = 200e3)
Im=489.86; //Moment of Inertia in mm4 (10mm diameter bar = 489.86)
Gm=Em;
Jm=2*Im;
Am=[78.4787]; //Area in mm2 (10mm diameter bar = 78.4787)

//If any properties change for each element, use matrix input

//Create Properties Matrix
[Prop]=Properties(Am,Em,Im,Gm,Jm,conec);
