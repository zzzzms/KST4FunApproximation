function x=myOMP(A,b,NIt)
% This is a  version of OMP based on Tropp written in Oct., 2007.  It
%was modified by Dr. Ming-Jun Lai on Sept. 2010.  
xsize=size(A,1); ysize=size(A,2);
y2=b;
IND=zeros(ysize,1); y=zeros(ysize,1);
for i=1:NIt
xa=abs(A'*b);
mx=max(xa)*0.995;
IND=IND+(xa>mx); IND=IND>0;
A1=A(:,IND); y=A1\b;
x=zeros(ysize,1); x(IND)=y;
    z=A*x; b=b-z;
   err=norm(b,inf);
   if err<1e-8 
       break
   end
end
%N=sum(IND);
 A1=A(:,IND); y=A1\y2;
 x=zeros(ysize,1); x(IND)=y; 
 norm(A*x-y2,inf)