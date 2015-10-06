// Nemo Frame preliminary design script
//
// Finite Element Aproach using 2D beam elements
//.First step to 3D fem script, just to get the hang of it
//
// Worth noting as a reference the Article
// "Overview of a method that can help select and refine the optimum human-powered vehicle frame design" - Don Chain
// "Human Power" Magazine, Spring 1989, Vol7 No.3 - p.6
//
// Diego Montero; Fernando Valentini; Gustavo Violato; 
// First release: Oct. 2015 
clc; xdel(winsid()); clear;


//----------------USER DEFINED VARIABLES----------------------//

//Material properties & Geometry
Em=200e3; //Young Modulus in MPa (1020 Steel = 200e3)
Im=489.86; //Moment of Inertia in mm4 (10mm diameter bar = 489.86)

 //Area in mm2 (10mm diameter bar = 78.4787)
Am=[78.4787; 
78.4787;
78.4787
];
//If any properties change for each element, use matrix input

//Nodes coordinates, each line is a node n located at x & y
coord = [   0 0
            100 0
            50 75
];

//Conectivity Matrix, each line is an element connecting nodes a & b
conec = [1 2
2 3
1 3
];

//Forces, applied to matching DOF
    F(8)=-30000; //in this case F in y of node 3 (DOF 8)

//Restricted DOFs (also for boundary conditions)
glrest=[1; 2; 5];

//----------------PROGRAM----------------------//

//Elements, Nodes & Degrees of Freedom
elem=size(conec,1);
nos=size(coord,1);
ngl = 3*nos;

//Global Matrices inicialization
Ul = zeros(4,1);
Re = zeros(4,1);
Fg = resize_matrix(F,ngl);
Kg = zeros(ngl,ngl);
N=zeros(1,elem);
V=zeros(1,elem);
M=zeros(2,elem);
PP = zeros(3,1);

//Local and Global Stiffness Matrices (2D case)
for e=1:elem
    //Calls for element material properties;
   if size(Em,1)>1
       E=Em(elem)
   else E=Em(1); end
   if size(Am,1)>1
       A=Am(elem)
   else A=Am(1); end
   if size(Im,1)>1
       I=Im(elem)
   else I=Im(1); end
     //Define coordinates for a given element
    n1 = conec(e,1);
    n2 = conec(e,2);
    x1 = coord(n1,1);
    y1 = coord(n1,2);
    x2 = coord(n2,1);
    y2 = coord(n2,2);
    //Euclidian Norm between nodes
    t = sqrt((x2-x1)^2+(y2-y1)^2);
    //Local Stiffness Matrix of the element (BEAM)
    K = [E*A/t 0 0 -(E*A/t) 0 0;
    0 12*E*I/(t^3) 6*E*I/(t^2) 0 -12*E*I/(t^3) 6*E*I/(t^2);
    0 6*E*I/(t^2) 4*E*I/(t) 0 -6*E*I/(t^2) 2*E*I/(t);
    -E*A/t 0 0 E*A/t 0 0;
    0 -12*E*I/(t^3) -6*E*I/(t^2) 0 12*E*I/(t^3) -6*E*I/(t^2);
    0 6*E*I/(t^2) 2*E*I/(t) 0 -6*E*I/(t^2) 4*E*I/(t)];
    //Stores local stiffness of each element
    Klocal(:,:,e)=K;
    //Stiffness Matrix rotation
    c = (x2-x1)/t;
    s = (y2-y1)/t;
    R = [ c s 0 0 0 0; -s c 0 0 0 0; 0 0 1 0 0 0;
        0 0 0 c s 0; 0 0 0 -s c 0; 0 0 0 0 0 1];
    //Stores local Rotation Matrix
    Rlocal(:,:,e)=R;
    //Rotates K to insert into Global Matrix
    K = R'* K * R;
    gl = [3*n1-2 3*n1-1 3*n1 3*n2-2 3*n2-1 3*n2];
    Kg(gl,gl) = Kg(gl,gl) + K;
end
Kgl = Kg; //Saves original Stiffness Matrix to Kgl (no DOF contraints)

//Limits the chosen DOFs in the Global Stiffness Matrix
nglr=size(glrest,1);
for i=1:nglr
    Kg(glrest(i),:) = 0;
    Kg(:,glrest(i)) = 0;
    Kg(glrest(i),glrest(i)) = 1;
end

//Nodes Displacement Matrix (mm)
Ug = Kg\Fg;
//Nodes Global reactions (N)
Re = Kgl*Ug;

//Local Calculation for Normal, Shear and Momentum
for e=1:elem
    //Calls Local R and K Matrices
    n1 = conec(e,1);
    n2 = conec(e,2);
    K=Klocal(:,:,e);
    R=Rlocal(:,:,e);
    if size(Em,1)>1
    E=Em(elem)
    else E=Em(1); end
    if size(Am,1)>1
    A=Am(elem)
    else A=Am(1); end
    if size(Im,1)>1
    I=Im(elem)
    else I=Im(1); end
    
    //Local Displacements from Global
    Ul(1,1) = Ug(3*n1-2);
    Ul(2,1) = Ug(3*n1-1);
    Ul(3,1) = Ug(3*n1);
    Ul(4,1) = Ug(3*n2-2);
    Ul(5,1) = Ug(3*n2-1);
    Ul(6,1) = Ug(3*n2);
    
    //Force=Stiffness*(Rotated Diplacements)
    P=K*R*Ul;  
   
   //Internal forces per node
    N(1,e)=P(1,1);
    N(2,e)=P(4,1);
    V(1,e)=P(2,1);
    V(2,e)=P(5,1);
    M(1,e)=P(3,1);
    M(2,e)=P(6,1);

//Calculate Stresses
Sa(:,e)=N(:,e)/A;
Sv(:,e)=V(:,e)/A;
Sm(:,e)=M(:,e)/I;

//Von-mises stress
Snn(n1,e)=sqrt((Sm(1,e)+Sa(1,e)).^2 + 3*Sv(1,e).^2)
Snn(n2,e)=sqrt((Sm(2,e)+Sa(2,e)).^2 + 3*Sv(2,e).^2)
end

//Organize Stress per node, always getting higher stresses
for n=1:nos
    [m,k]=max(abs(Snn(n,:)))
    S(n)=Snn(n,k);
    Udef(n,:) = [Ug(3*n-2) Ug(3*n-1)];
    Ulin(n)=sqrt(Ug(3*n-2)^2+Ug(3*n-1)^2);
    Utor(n)=180*Ug(3*n)/%pi;
    Results(n,1)=n;
end
//Organize and display Results;
Results(:,2)=S;
Results(:,3)=Ulin;
Results(:,4)=Utor;
printf("  Node   Stress(MPa) Lin.Disp.(mm) Ang.Disp.(°)")
disp(Results);

//Plot truss

function plottruss(T, P, _color)
//---------------------------------------------------
// 
//  T is the connectivity between nodes (or conec)
//  P is the nodes coordinates (or coord)

//To simplify:
NodeLabels=0;
EltLabels=0;
//Though these variables might be in the function if needed.

Net = size(T,1);
_3D_problem = (size(P,2)==3);

if _3D_problem then
  for ie = 1:Net
    XY = P(T(ie,:),:);
    X = [[XY(:,1)]; XY(1,1)];
    Y = [[XY(:,2)]; XY(1,2)];
    Z = [[XY(:,3)]; XY(1,3)];
    param3d1(X,Y,list(Z,color(_color)));  
    if (EltLabels)
      x = mean(XY(:,1));
      y = mean(XY(:,2));
      z = mean(XY(:,3));
      xstring(x,y,z,string(i));
    end
  end 

  if (NodeLabels) then
    Np = size(P,1);
    for i=1:Np
      xstring(P(i,1),P(i,2),P(i,3),string(i));
    end
  end

  delta_x = max(P(:,1))-min(P(:,1));
  delta_y = max(P(:,2))-min(P(:,2));
  delta_z = max(P(:,3))-min(P(:,3));
  x_max = max(P(:,1));
  x_min = min(P(:,1));
  y_max = max(P(:,2));
  y_min = min(P(:,2));
  z_max = max(P(:,3));
  z_min = min(P(:,3));

  f = gca();
  f.data_bounds = [x_min-0.1*delta_x, y_min-0.1*delta_y, z_min-0.1*delta_z;
                   x_max+0.1*delta_x, y_max+0.1*delta_y, z_max+0.1*delta_z];
else
  for ie = 1:Net
    XY = P(T(ie,:),:);
    X = [[XY(:,1)]; XY(1,1)];
    Y = [[XY(:,2)]; XY(1,2)];
    plot(X,Y,_color);  
    if (EltLabels)
      x = mean(XY(:,1));
      y = mean(XY(:,2));
      xstring(x,y,string(ie));
    end
  end 

  if (NodeLabels) then
    Np = size(P,1);
    for i=1:Np
      xstring(P(i,1),P(i,2),string(i));
    end
  end

  delta_x = max(P(:,1))-min(P(:,1));
  delta_y = max(P(:,2))-min(P(:,2));
  x_max = max(P(:,1));
  x_min = min(P(:,1));
  y_max = max(P(:,2));
  y_min = min(P(:,2));

  f = gca();
  f.data_bounds = [x_min-0.1*delta_x, y_min-0.1*delta_y; x_max+0.1*delta_x, y_max+0.1*delta_y];
end
endfunction
plottruss(conec,coord,'green');
coord=coord+10*Udef;
plottruss(conec,coord,'red');
legends(['before loading','after loading'], [color('green') color('red')], 'ur');

// Change log

// 2015-10-04 - Diego Montero
// ----------------------------
// First release, working for 2D truss with simple plot in a single file.
