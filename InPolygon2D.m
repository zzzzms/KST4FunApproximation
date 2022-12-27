function [ in, on ] = InPolygon2D( X, Y, XV, YV, tol )
%INPOLYGON2D Summary of this function goes here
%   Detailed explanation goes here
% Written by Yidong Xu.
% Input format follows that of Matlab's built-in function inpolygon.

x_min=min(XV);
x_max=max(XV);
y_min=min(YV);
y_max=max(YV);

% E: the coordinates of the two endpoints.
E=zeros(length(XV),4);
% AnalyzeEdges: 
% The first 4 components are x_min, x_max, y_min, and y_max.
% The 5th component is sqrt((x2-x1)^2+(y2-y1)^2).
AnalyzeEdges=zeros(length(XV),5);
num_of_edges=0;

for i=1:length(XV)-1
    if isnan(XV(i)) || isnan(XV(i+1))
       continue; 
    end
    num_of_edges=num_of_edges+1;
    E(num_of_edges,:)=[XV(i),YV(i),XV(i+1),YV(i+1)];
    AnalyzeEdges(num_of_edges,:)=[min(XV(i),XV(i+1)),max(XV(i),XV(i+1)),...
        min(YV(i),YV(i+1)),max(YV(i),YV(i+1)),...
        sqrt((XV(i+1)-XV(i))^2+(YV(i+1)-YV(i))^2)];
end

E(num_of_edges+1:end,:)=[];
AnalyzeEdges(num_of_edges+1:end,:)=[];
% index=1:num_of_edges;
% E=E(index,:);
% AnalyzeEdges=AnalyzeEdges(index,:);

on=false(length(X),1);
in=on;

FindIndex=find(X>=x_min-tol & X<=x_max+tol & Y>=y_min-tol & Y<=y_max+tol);
n=length(FindIndex);

if n==0
    return;
end

% block_size for vectorization.
block_size=1000;
num_of_testing_points_in_each_block=floor(block_size/num_of_edges);
if num_of_testing_points_in_each_block==0
    num_of_testing_points_in_each_block=1;
end
block_num=ceil(n/num_of_testing_points_in_each_block);
if block_num==1
    num_of_testing_points_in_each_block=n;
end

X1=repmat(E(:,1),num_of_testing_points_in_each_block,1);
Y1=repmat(E(:,2),num_of_testing_points_in_each_block,1);
X2=repmat(E(:,3),num_of_testing_points_in_each_block,1);
Y2=repmat(E(:,4),num_of_testing_points_in_each_block,1);
AnalyzeEdges1=repmat(AnalyzeEdges(:,1),num_of_testing_points_in_each_block,1);
AnalyzeEdges2=repmat(AnalyzeEdges(:,2),num_of_testing_points_in_each_block,1);
AnalyzeEdges3=repmat(AnalyzeEdges(:,3),num_of_testing_points_in_each_block,1);
AnalyzeEdges4=repmat(AnalyzeEdges(:,4),num_of_testing_points_in_each_block,1);
AnalyzeEdges5=repmat(AnalyzeEdges(:,5),num_of_testing_points_in_each_block,1);


for i=1:block_num-1
    start=(i-1)*num_of_testing_points_in_each_block;
    index=start+1:start+num_of_testing_points_in_each_block;
    x=reshape(repmat(X(FindIndex(index))',num_of_edges,1),[],1);
    y=reshape(repmat(Y(FindIndex(index))',num_of_edges,1),[],1);
    
    status=AnalyzeEdges1<=x+tol & AnalyzeEdges2>=x-tol & AnalyzeEdges3<=y+tol ...
        & AnalyzeEdges4>=y-tol & abs((X2-X1).*(y-Y1)-(Y2-Y1).*(x-X1))<=tol*AnalyzeEdges5;
    
    on(FindIndex(start+1:start+num_of_testing_points_in_each_block))=max(reshape(status,[],num_of_testing_points_in_each_block));
    
end

start=(block_num-1)*num_of_testing_points_in_each_block;
num_of_points_remaining=n-start;
index=start+1:n;
x=reshape(repmat(X(FindIndex(index))',num_of_edges,1),[],1);
y=reshape(repmat(Y(FindIndex(index))',num_of_edges,1),[],1);
index=1:num_of_points_remaining*num_of_edges;
status=AnalyzeEdges1(index)<=x+tol & AnalyzeEdges2(index)>=x-tol ...
    & AnalyzeEdges3(index)<=y+tol & AnalyzeEdges4(index)>=y-tol ...
    & abs((X2(index)-X1(index)).*(y-Y1(index))-(Y2(index)-Y1(index)).*(x-X1(index)))<=tol*AnalyzeEdges5(index);

on(FindIndex(start+1:n))=max(reshape(status,[],num_of_points_remaining));


in(FindIndex)=on(FindIndex) | inpolygon(X(FindIndex),Y(FindIndex),XV,YV);

return;












% Old code below. Don't delete.

x_min=min(XV);
x_max=max(XV);
y_min=min(YV);
y_max=max(YV);

% E: the coordinates of the two endpoints.
E=zeros(length(XV),4);
% AnalyzeEdges: 
% The first 4 components are x_min, x_max, y_min, and y_max.
% The 5th component is sqrt((x2-x1)^2+(y2-y1)^2).
AnalyzeEdges=zeros(length(XV),5);
num_of_edges=0;

for i=1:length(XV)-1
    if isnan(XV(i)) || isnan(XV(i+1))
       continue; 
    end
    num_of_edges=num_of_edges+1;
    E(num_of_edges,:)=[XV(i),YV(i),XV(i+1),YV(i+1)];
    AnalyzeEdges(num_of_edges,:)=[min(XV(i),XV(i+1)),max(XV(i),XV(i+1)),...
        min(YV(i),YV(i+1)),max(YV(i),YV(i+1)),...
        sqrt((XV(i+1)-XV(i))^2+(YV(i+1)-YV(i))^2)];
end

E(num_of_edges+1:end,:)=[];
AnalyzeEdges(num_of_edges+1:end,:)=[];
% index=1:num_of_edges;
% E=E(index,:);
% AnalyzeEdges=AnalyzeEdges(index,:);

% n=length(X);
on=false(length(X),1);

% FindIndex1=find(X>=x_min-tol & X<=x_max+tol & Y>=y_min-tol & Y<=y_max+tol);
% 
% for i=1:num_of_edges
%     x1=E(i,1);
%     y1=E(i,2);
%     x2=E(i,3);
%     y2=E(i,4);
%     
%     FindIndex2=X(FindIndex1)>=min(x1,x2)-tol & X(FindIndex1)<=max(x1,x2)+tol ...
%         & Y(FindIndex1)>=min(y1,y2)-tol & Y(FindIndex1)<=max(y1,y2)+tol;
%     FindIndex2=FindIndex1(FindIndex2);
%     FindIndex3=abs((x2-x1).*(Y(FindIndex2)-y1)-(y2-y1).*(X(FindIndex2)-x1))<=tol;
%     FindIndex3=FindIndex2(FindIndex3);
%     on(FindIndex3)=true;
% end

FindIndex1=find(X>=x_min-tol & X<=x_max+tol & Y>=y_min-tol & Y<=y_max+tol);

for i=1:length(FindIndex1)
    x=X(FindIndex1(i));
    y=Y(FindIndex1(i));
%     if x<x_min-tol || x>x_max+tol || y<y_min-tol || y>y_max+tol
%         continue;
%     end
    
    FindIndex2=AnalyzeEdges(:,1)<=x+tol & AnalyzeEdges(:,2)>=x-tol & AnalyzeEdges(:,3)<=y+tol & AnalyzeEdges(:,4)>=y-tol;
    FindIndex2=find(abs((E(FindIndex2,3)-E(FindIndex2,1)).*(y-E(FindIndex2,2)) ...
        -(E(FindIndex2,4)-E(FindIndex2,2)).*(x-E(FindIndex2,1)))<=tol*AnalyzeEdges(FindIndex2,5),1);
    if ~isempty(FindIndex2)
        on(FindIndex1(i))=true;
    end
    
end

in=on | inpolygon(X,Y,XV,YV);


end

