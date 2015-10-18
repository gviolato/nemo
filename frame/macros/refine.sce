function [Ref_conec,Ref_coord,Ref_Prop]=refine(conec,coord,Prop,Ref)

//This function refines the Mesh (Coordinates and Conectivity)to a Ref value,
//Also expands the Properties matrix.

Ref_coord=coord;
Ref_conec=conec;
Ref_Prop=Prop;

//Do nothing if Ref = 0
if Ref~=0 then
[elem,sc]=size(conec);
[nos,sc]=size(coord);

//Brings a R2 problem to R3
if sc==2 then
   Ref_coord=resize_matrix(Ref_coord,nos,3);
end

for e=1:elem
    clear NewCoord;
    n1 = Ref_conec(e,1);
    n2 = Ref_conec(e,2);
    x1 = Ref_coord(n1,1);
    y1 = Ref_coord(n1,2);
    z1 = Ref_coord(n1,3);
    x2 = Ref_coord(n2,1);
    y2 = Ref_coord(n2,2);
    z2 = Ref_coord(n2,3);
    xin=x2-x1;
    yin=y2-y1;
    zin=z2-z1;
    
    //Calculate the element size to define the number of sub-sections
    t = sqrt((xin)^2+(yin)^2+(zin)^2);
    new=round(t/Ref);
    xadd=xin/new;
    yadd=yin/new;
    zadd=zin/new;
    
    //Adds the subsections to NewCoord
    NewCoord(1,:)=[x1+xadd y1+yadd z1+zadd];
    for r=2:(new-1)
        NewCoord(r,:)=NewCoord(r-1,:)+[xadd yadd zadd];
    end
    
    //Inputs NewCoord into coordinates matrix
    [a1,b1]=size(Ref_coord);
    Ref_coord((a1+1):(a1+new-1),:)=NewCoord(:,:);
    
    //Corrects conectivity matrix
    Ref_conec(e,:)=[n1 (a1+1)];
    for r=2:(new)
        [a2,b2]=size(Ref_conec);
        Ref_conec(a2+1,:)=[a1+(r-1) (a1+r)];
    //And material properties
        Ref_Prop(a2+1,:)=Prop(e,:);
    end
    Ref_conec(a2+1,2)=n2;
end

//Goes back to R2 if needed
if sc==2 then
   [a3,b3]=size(Ref_coord);
   Ref_coord=resize_matrix(Ref_coord,a3,(sc));
end
end

endfunction
