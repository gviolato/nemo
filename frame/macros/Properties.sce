function [Prop]=Properties(Am,Em,Im,Gm,Jm,conec)

//This function compiles the material properties in a matrix,
//main purpose is to simplify the problem inputs

[elem,bb]=size(conec);
//Correct Material Properties Matrices / warn if wrong input was given
warn=0;
   if size(Am,1)==1
       A(1:elem)=Am;
   elseif size(Am,1)>1
   warn=1; end
   if size(Em,1)==1
       E(1:elem)=Em;
   elseif size(Em,1)>1
   warn=1; end
      if size(Im,1)==1
       I(1:elem)=Im;
   elseif size(Im,1)>1
   warn=1; end
         if size(Gm,1)==1
       G(1:elem)=Gm;
   elseif size(Gm,1)>1
   warn=1; end
         if size(Jm,1)==1
       J(1:elem)=Jm;
   elseif size(Jm,1)>1
   warn=1; end
   if warn==1 then
       printf("WARNING: BAD MATERIAL PROPERTIES INPUT")
   end

Prop=[A E I G J];

   endfunction
