%SPLINEFITTING3D 
% Written by Yidong Xu under Dr. Ming-Jun Lai's supervision, 2014--2019.
%Please acknowledge Dr. Ming-Jun Lai
%and Yidong Xu for their contribution and generocity.
global Factorials;
n=(d+2)*(d+1)/2;m=size(T,1);
TR=triangulation(T,V);
[TriangleIndex,Bary]=pointLocation(TR,X,Y);
FindIndex=find(isnan(TriangleIndex));
if ~isempty(FindIndex)
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

eps1=m^(2);
[row,column]=size(H); L=A'*A+lambda*K;
%c=ConstrainedLeastSquareKKT(A,W,K,H,lambda,tol);
U=H'*H;
dA = decomposition(U+1/eps1*L,'auto');
