function [Results,Udef]=solve_truss(conec, coord, Prop, F)
    
//Main program, solves a 2D or 3D truss with beam elements.
//Even though there is a single function, the function reads the
//problem dimention and execute the accordinly routine.
//Outputs are von-mises stresses and linear displacementes

//Elements, Nodes & dimention
[nos,sc]=size(coord);
elem=size(conec,1);

//Loads Properties
A=Prop(:,1); E=Prop(:,2); I=Prop(:,3); G=Prop(:,4); J=Prop(:,5);

//Assures that restricted DOF were entered correctly
[dbg,nglr]=size(glrest) 
if nglr==1 
    glrest=glrest';
end

//----------------PROGRAM 3D----------------------//
if sc==3 then
//Degrees of Freedom
ngl = 6*nos;

//Global Matrices inicialization
Ul = zeros(6,1);
Re = zeros(6,1);
Fg = resize_matrix(F,ngl);
Kg = zeros(ngl,ngl);
N=zeros(2,elem);
V=zeros(2,elem);
M=zeros(2,elem);

//Local and Global Stiffness Matrices (2D case)
for e=1:elem
    //Define coordinates for a given element
    n1 = conec(e,1);
    n2 = conec(e,2);
    x1 = coord(n1,1);
    y1 = coord(n1,2);
    z1 = coord(n1,3);
    x2 = coord(n2,1);
    y2 = coord(n2,2);
    z2 = coord(n2,3);
    //Euclidian Norm between nodes
    t = sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
    a1=E(e)*A(e)/t;
    a2=12*E(e)*I(e)/(t^3);
    a3=6*E(e)*I(e)/(t^2);
    a4=4*E(e)*I(e)/t;
    a5=G(e)*J(e)/t;
    a6=a4/2;
    //Local Stiffness Matrix of the element (3DBEAM with Iy=Iz)
    K = [a1 0 0 0 0 0     -a1 0 0 0 0 0;
         0 a2 0 0 0 a3    0 -a2 0 0 0 a3;
         0 0 a2 0 -a3 0   0 0 -a2 0 -a3 0;
         0 0 0 a5 0 0     0 0 0 -a5 0 0;
         0 0 -a3 0 a4 0   0 0 a3 0 a6 0;
         0 a3 0 0 0 a4    0 -a3 0 0 0 a6;
                                       
         -a1 0 0 0 0 0    a1 0 0 0 0 0;
         0 -a2 0 0 0 -a3  0 a2 0 0 0 -a3;
         0 0 a2 0 a3 0    0 0 a2 0 a3 0;
         0 0 0 -a5 0 0    0 0 0 a5 0 0;
         0 0 -a3 0 a6 0   0 0 a3 0 a4 0;
         0 a3 0 0 0 a6    0 -a3 0 0 0 a4];
    //Stores local stiffness of each element
    Klocal(:,:,e)=K;
    
    //Stiffness Matrix rotation
    ll = (x2-x1)/t;
    mm = (y2-y1)/t;
    nn = (z2-z1)/t;
    At=[ll mm nn];
    Bt=cross([0 0 1],At);
    NBt=sqrt(sum(Bt.^2));
    Bt=Bt/NBt;
    Ct=cross(At,Bt);
    Tag(1,1:3)=At;
    Tag(2,1:3)=Bt;
    Tag(3,1:3)=Ct;
    Tla=[1 0 0;
         0 ll mm;
         0 -mm ll];
    LB=Tla*Tag;
    R=zeros(12,12);
    R(1:3,1:3)=LB;
    R(4:6,4:6)=LB;
    R(7:9,7:9)=LB;
    R(10:12,10:12)=LB;
    
    //Stores local Rotation Matrix
    Rlocal(:,:,e)=R;
    //Rotates K to insert into Global Matrix
    K = R'* K * R;
    gl = [6*n1-5 6*n1-4 6*n1-3 6*n1-2 6*n1-1 6*n1     6*n2-5 6*n2-4 6*n2-3 6*n2-2 6*n2-1 6*n2];
    Kg(gl,gl) = Kg(gl,gl) + K;
end
Kgl = Kg; //Saves original Stiffness Matrix to Kgl (no DOF contraints)

//Limits the chosen DOFs in the Global Stiffness Matrix
[dbg,nglr]=size(glrest);
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
    K2=Klocal(:,:,e);
    R=Rlocal(:,:,e);
    //Local Displacements from Global
    Ul(1,1) = Ug(6*n1-5);
    Ul(2,1) = Ug(6*n1-4);
    Ul(3,1) = Ug(6*n1-3);
    Ul(4,1) = Ug(6*n1-2);
    Ul(5,1) = Ug(6*n1-1);
    Ul(6,1) = Ug(6*n1);
    Ul(7,1) = Ug(6*n2-5);
    Ul(8,1) = Ug(6*n2-4);
    Ul(9,1) = Ug(6*n2-3);
    Ul(10,1) = Ug(6*n2-2);
    Ul(11,1) = Ug(6*n2-1);
    Ul(12,1) = Ug(6*n2);
    
    //Force=Stiffness*(Rotated Diplacements)
    P=K2*R*Ul;  
   
   //Internal forces per node
    N(1,e)=P(1,1);
    N(2,e)=P(7,1);
    V1(1,e)=P(2,1);
    V1(2,e)=P(8,1);
    V2(1,e)=P(3,1);
    V2(2,e)=P(9,1);
    T(1,e)=P(4,1);
    T(2,e)=P(10,1);
    M1(1,e)=P(5,1);
    M1(2,e)=P(11,1);
    M2(1,e)=P(6,1);
    M2(2,e)=P(12,1);

V=(V1.^2+V2.^2).^(1/2);
M=(M1.^2+M2.^2).^(1/2);

//Calculate Stresses
Sa(:,e)=N(:,e)/A(e);
Sv(:,e)=V(:,e)/A(e);
Sm(:,e)=M(:,e)/I(e);
St(:,e)=T(:,e)/J(e);

//Von-mises stress
Snn(n1,e)=sqrt((Sm(1,e)+Sa(1,e)).^2 + 3*(Sv(1,e)+St(1,e)).^2)
Snn(n2,e)=sqrt((Sm(2,e)+Sa(2,e)).^2 + 3*(Sv(2,e)+St(2,e)).^2)
end

//Organize Stress per node, always getting higher stresses
for n=1:nos
    [aaa,k]=max(abs(Snn(n,:)))
    S(n)=Snn(n,k);
    Udef(n,:) = [Ug(6*n-5) Ug(6*n-4) Ug(6*n-3)];
    Ulin(n)=sqrt(sum(Udef(n,:).^2));
    Utor(n)=180*Ug(3*n)/%pi;
    Results(n,1)=n;
end
//Organize and display Results;
Results(:,2)=S;
Results(:,3)=Ulin;
end

    //----------------PROGRAM 2D----------------------//

if sc==2 then

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
     //Define coordinates for a given element
    n1 = conec(e,1);
    n2 = conec(e,2);
    x1 = coord(n1,1);
    y1 = coord(n1,2);
    x2 = coord(n2,1);
    y2 = coord(n2,2);
    //Euclidian Norm between nodes
    t = sqrt((x2-x1)^2+(y2-y1)^2);
    a1=E(e)*A(e)/t;
    a2=12*E(e)*I(e)/(t^3);
    a3=6*E(e)*I(e)/(t^2);
    a4=4*E(e)*I(e)/t;
    a6=2*E(e)*I(e)/t;
    //Local Stiffness Matrix of the element (BEAM)
    K = [a1 0 0    -a1 0 0;
         0 a2 a3    0 -a2 a3;
         0 a3 a4    0 -a3 a6;
                            
         -a1 0 0    a1 0 0;
         0 -a2 -a3  0 a2 -a3;
         0 a3 a6    0 -a3 a4];
    //Stores local stiffness of each element
    Klocal(:,:,e)=K;
    //Stiffness Matrix rotation
    c = (x2-x1)/t;
    s = (y2-y1)/t;
    R = [c s 0    0 0 0;
         -s c 0   0 0 0;
         0 0 1    0 0 0;
                        
         0 0 0    c s 0;
         0 0 0    -s c 0;
         0 0 0    0 0 1];
    //Stores local Rotation Matrix
    Rlocal(:,:,e)=R;
    //Rotates K to insert into Global Matrix
    K = R'* K * R;
    gl = [3*n1-2 3*n1-1 3*n1 3*n2-2 3*n2-1 3*n2];
    Kg(gl,gl) = Kg(gl,gl) + K;
end
Kgl = Kg; //Saves original Stiffness Matrix to Kgl (no DOF contraints)

//Limits the chosen DOFs in the Global Stiffness Matrix
[dbg,nglr]=size(glrest);
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
Sa(:,e)=N(:,e)/A(e);
Sv(:,e)=V(:,e)/A(e);
Sm(:,e)=M(:,e)/I(e);

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
//    Utor(n)=180*Ug(3*n)/%pi;
    Results(n,1)=n;
end
//Organize and display Results;
Results(:,2)=S;
Results(:,3)=Ulin;
//Results(:,4)=Utor;
end
endfunction
