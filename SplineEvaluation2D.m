function Val= SplineEvaluation2D(V,T,Analyze,c,d,X,Y,tol)
% This code is written by Yidong Xu under Dr. Ming-Jun Lai's supervision 
% during 2014-2019. 
%Please acknowledge Dr. Ming-Jun Lai
%and Yidong Xu for their contribution and generocity.
% Only one polynomial and many points
% If the i-th point is not in any tetrahedron, Val(i)=NaN.

Val=nan(length(X),1);

TR=triangulation(T,V);
[TriangleIndex,Bary]=pointLocation(TR,X,Y);
FindIndex=find(isnan(TriangleIndex));

if ~isempty(FindIndex)
    [TriangleIndex(FindIndex),Bary(FindIndex,:)]=LocatePoints2D(V,T,Analyze,X(FindIndex),Y(FindIndex),tol);   
end

FindIndex=find(~isnan(TriangleIndex));
if isempty(FindIndex)
    return;
end

if d==0
%    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    Val=c(TriangleIndex);
    return;
end

p=length(FindIndex);
Bary=Bary';


global C1Index;
global C2Index;
global C3Index;

% n=(d+2)*(d+1)/2;
% Val_temp=reshape(c,n,[]);
m=size(T,1);
Val_temp=reshape(c,[],m);
%Val=reshape(Val(:,TriangleIndex),[],1);
Val_temp=Val_temp(:,TriangleIndex(FindIndex));

for i=d:-1:1
%     index=1:(i+2)*(i+1)/2;
    c1=reshape(Val_temp(repmat(C1Index{i,1},p,1)),[],p);
    c2=reshape(Val_temp(repmat(C2Index{i,1},p,1)),[],p);
    c3=reshape(Val_temp(repmat(C3Index{i,1},p,1)),[],p);
    Val_temp=c1*spdiags(Bary(1,FindIndex)',0,p,p)+c2*spdiags(Bary(2,FindIndex)',0,p,p)...
        +c3*spdiags(Bary(3,FindIndex)',0,p,p);   
end
Val(FindIndex)=Val_temp;
end

