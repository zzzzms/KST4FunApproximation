function c= SplineFitting2D(V,T,Analyze,d,r, X, Y, Z,lambda,tol,opts )
% This matlab code is written by Yidong Xu under Dr. Ming-Jun Lai's
% supervision in 2019.  It is rewritten by Dr. Ming-Jun Lai on May 18,
% 2022.
% opts: struct with fields:
% NIterations,x_partition_num,y_partition_num,Boundary,Anticycle,InfNum

if isfield(opts,'NIterations')
    NIterations=opts.NIterations;
end
if isfield(opts,'x_partition_num')
    x_partition_num=opts.x_partition_num;
end
if isfield(opts,'y_partition_num')
    y_partition_num=opts.y_partition_num;
end
if isfield(opts,'Boundary')
    Boundary=opts.Boundary;
end
if isfield(opts,'Anticycle')
    Anticycle=opts.Anticycle;
end
if isfield(opts,'InfNum')
    InfNum=opts.InfNum;
end
% if isfield(opts,'MaxDimension')
%     MaxDimension=opts.MaxDimension;
% end

global Factorials;

n=(d+2)*(d+1)/2;m=size(T,1);
TR=triangulation(T,V);
[TriangleIndex,Bary]=pointLocation(TR,X,Y);
FindIndex=find(isnan(TriangleIndex));
% FindIndex=find(isnan(TriangleIndex));
% length(FindIndex);
% [X(FindIndex),Y(FindIndex)]
if ~isempty(FindIndex)
%     disp('***************************');
%     length(FindIndex)
[TriangleIndex(FindIndex),Bary(FindIndex,:)]=LocatePoints2D(V,T,Analyze,X(FindIndex),Y(FindIndex),tol);
    
FindIndex2=FindIndex(isnan(TriangleIndex(FindIndex)));
    if ~isempty(FindIndex2)
        X(FindIndex2)=[];
        Y(FindIndex2)=[];
        Z(FindIndex2)=[];
        TriangleIndex(FindIndex2)=[];
        Bary(FindIndex2,:)=[];
    end  
end

p=length(X);

% Note that barypowers(1,1,1)=b1^0=1, barypowers(1,1,2)=b1^1=b1, and so on.
barypowers=ones(p,3);
if d>=1
    barypowers(:,:,2)=Bary;
end
for i=3:d+1
    barypowers(:,:,i)=barypowers(:,:,i-1).*Bary;
end

BCoef=zeros(p,n);
index=0;
for i=d:-1:0
   for j=d-i:-1:0
       index=index+1; k=d-i-j;
   BCoef(:,index)=(Factorials(d+1)/Factorials(i+1)/Factorials(j+1)/Factorials(k+1))...
            *barypowers(:,1,i+1).*barypowers(:,2,j+1).*barypowers(:,3,k+1);
    end
end
A=sparse(repmat((1:p)',1,n),repmat((TriangleIndex-1)*n,1,n)+repmat(1:n,p,1),BCoef,p,n*m);
H=Smoothness2D(V,T,d,r); % Smoothness matrix H
K=Energy2D(Analyze,d); % Energy matrix K
warning('off','MATLAB:rankDeficientMatrix');
  c=ConstrainedLeastSquareKKT(A,Z,K,H,lambda,tol);
warning('on','MATLAB:rankDeficientMatrix');
end


function x = ConstrainedLeastSquareKKT(A,b,K,H,lambda,tol)
b=sparse(b); [row,column]=size(H);
if isempty(H)
%     if row1==column1
%         x=[A;sparse(1,column1)]\[b;0];
%     else
%         x=A\b;
%     end
    x=[A'*A+lambda*K;sparse(1,column)]\[A'*b;0];
    return;
end
M=[A'*A+lambda*K,H';H,sparse(row,row);sparse(1,column+row)];
N=[A'*b;sparse(row+1,1)];
x=M\N;
if norm(M*x-N,inf)>tol
    disp('The solution may not be correct.');
end
x=x(1:column);
end


