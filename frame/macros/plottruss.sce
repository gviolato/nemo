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
    //Switching axes!
    ss=P(:,2);
P(:,2)=P(:,3);
P(:,3)=ss;
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
