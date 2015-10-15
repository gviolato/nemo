//3D Tretrahedron Problem

//Nodes coordinates, each line is a node n located at x, y & z
coord = [   0 0 0
            100 0 0
            100 100 0
            0 100 0
            50 50 70.71068
];

//Conectivity Matrix, each line is an element connecting nodes a & b
conec = [1 2
        2 3
        3 4
        4 1
        1 5
        2 5
        3 5
        4 5
];

//Forces, applied to matching DOF
    F(25)=-50000; //in this case F in x of node 5 (DOF 25)

//Restricted DOFs
glrest=[1; 2; 3; 6+2; 12+2; 18+2];

//Material properties & Geometry
Em=200e3; //Young Modulus in MPa (1020 Steel = 200e3)
Im=489.86; //Moment of Inertia in mm4 (10mm diameter bar = 489.86)
Gm=Em;
Jm=2*Im;
Am=[78.4787]; //Area in mm2 (10mm diameter bar = 78.4787)

//If any properties change for each element, use matrix input

//Create Properties Matrix
[Prop]=Properties(Am,Em,Im,Gm,Jm,conec);


