function [ T, AnalyzeTriangulation ] = AnalyzeTriangulation2D( V, T )
%ANALYZETRIANGULATION2D Summary of this function goes here
%   Detailed explanation goes here
% Written by Yidong Xu.

% Important: signed area
X=[V(T(:,1),1),V(T(:,2),1),V(T(:,3),1)];
Y=[V(T(:,1),2),V(T(:,2),2),V(T(:,3),2)];

% m=size(T,1);
% three_m=3*m;
% M=sparse([reshape(repmat(1:3:three_m,3,1),1,[]),repmat(2:3:three_m,1,3),repmat(3:3:three_m,1,3)],             [1:three_m,repmat([1:3:three_m,2:3:three_m,3:3:three_m],1,2)],           [ones(three_m,1);V(T(:,1),1);V(T(:,2),1);V(T(:,3),1);V(T(:,1),2);V(T(:,2),2);V(T(:,3),2)]);

area_determinant=(X(:,2)-X(:,1)).*(Y(:,3)-Y(:,1))-(X(:,3)-X(:,1)).*(Y(:,2)-Y(:,1));
% Make the vertices of triangles counterclockwise.
FindIndex=area_determinant<0;
temp=T(FindIndex,1);
T(FindIndex,1)=T(FindIndex,2);
T(FindIndex,2)=temp;
area_determinant(FindIndex)=abs(area_determinant(FindIndex));
% Update X and Y.
temp=X(FindIndex,1);
X(FindIndex,1)=X(FindIndex,2);
X(FindIndex,2)=temp;
temp=Y(FindIndex,1);
Y(FindIndex,1)=Y(FindIndex,2);
Y(FindIndex,2)=temp;

% AnalyzeTriangulation=[min(X,[],2),max(X,[],2),min(Y,[],2),max(Y,[],2),...
%     ((X(:,2)-X(:,1)).*(Y(:,3)-Y(:,1))-(X(:,3)-X(:,1)).*(Y(:,2)-Y(:,1)))/2,...
%     Y(:,2)-Y(:,3), Y(:,3)-Y(:,1), Y(:,1)-Y(:,2), X(:,3)-X(:,2), X(:,1)-X(:,3), X(:,2)-X(:,1)];
AnalyzeTriangulation=[min(X,[],2),max(X,[],2),min(Y,[],2),max(Y,[],2),area_determinant,...
    (Y(:,2)-Y(:,3))./area_determinant, (Y(:,3)-Y(:,1))./area_determinant, (Y(:,1)-Y(:,2))./area_determinant,...
    (X(:,3)-X(:,2))./area_determinant, (X(:,1)-X(:,3))./area_determinant, (X(:,2)-X(:,1))./area_determinant];

AnalyzeTriangulation(:,5)=AnalyzeTriangulation(:,5)/2;

% AnalyzeTriangulation(:,6)=AnalyzeTriangulation(:,6)./AnalyzeTriangulation(:,5)/2;
% AnalyzeTriangulation(:,7)=AnalyzeTriangulation(:,7)./AnalyzeTriangulation(:,5)/2;
% AnalyzeTriangulation(:,8)=AnalyzeTriangulation(:,8)./AnalyzeTriangulation(:,5)/2;
% AnalyzeTriangulation(:,9)=AnalyzeTriangulation(:,9)./AnalyzeTriangulation(:,5)/2;
% AnalyzeTriangulation(:,10)=AnalyzeTriangulation(:,10)./AnalyzeTriangulation(:,5)/2;
% AnalyzeTriangulation(:,11)=AnalyzeTriangulation(:,11)./AnalyzeTriangulation(:,5)/2;

% make the vertices of triangles counterclockwise
% FindIndex=find(AnalyzeTriangulation(:,5)<0);
% temp=T(FindIndex,1);
% T(FindIndex,1)=T(FindIndex,2);
% T(FindIndex,2)=temp;
% AnalyzeTriangulation(FindIndex,5)=abs(AnalyzeTriangulation(FindIndex,5));

% n=size(T,1);
% AnalyzeTriangulation=zeros(n,4);
% 
% for i=1:n
%     X=[V(T(i,1),1),V(T(i,2),1),V(T(i,3),1)];
%     Y=[V(T(i,1),2),V(T(i,2),2),V(T(i,3),2)];
%     AnalyzeTriangulation(i,:)=[min(X),max(X),min(Y),max(Y)];
%     
% end


end

