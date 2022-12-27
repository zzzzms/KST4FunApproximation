function [ TriangleIndex, Bary ] = LocatePoints2D( V, T, Analyze, X, Y, tol )
%LOCATEPOINTS2D Summary of this function goes here
%   Detailed explanation goes here
% Written by Yidong Xu.
% If TriangleIndex(i)=NaN, the i-th point is not in any triangle.
% Then the corresponding Bary(i,:)=[NaN,NaN,NaN].

n=length(X);
%m=size(T,1);
TriangleIndex=nan(n,1);
Bary=nan(n,3);

for i=1:n
    x=X(i);
    y=Y(i);
    index=find(Analyze(:,1)<=x+tol & Analyze(:,2)>=x-tol & Analyze(:,3)<=y+tol & Analyze(:,4)>=y-tol);
    if isempty(index)
        continue;
    end
    t1=T(index,1);
    t2=T(index,2);
    t3=T(index,3);
    x1=V(t1,1);
    y1=V(t1,2);
    x2=V(t2,1);
    y2=V(t2,2);
    x3=V(t3,1);
    y3=V(t3,2);
    % the three sides of the triangle
%     FindIndex=((x2-x1).*(y3-y1)-(y2-y1).*(x3-x1)).*((x2-x1).*(Y(i)-y1)-(y2-y1).*(X(i)-x1))>=-tol & ((x3-x2).*(y1-y2)-(y3-y2).*(x1-x2)).*((x3-x2).*(Y(i)-y2)-(y3-y2).*(X(i)-x2))>=-tol & ((x3-x1).*(y2-y1)-(y3-y1).*(x2-x1)).*((x3-x1).*(Y(i)-y1)-(y3-y1).*(X(i)-x1))>=-tol;
    
    % the signed area
    FindIndex=find((x2-x1).*(y-y1)-(x-x1).*(y2-y1)>=-tol*sqrt((x2-x1).^2+(y2-y1).^2) ...
        & (x3-x2).*(y-y2)-(x-x2).*(y3-y2)>=-tol*sqrt((x3-x2).^2+(y3-y2).^2) ...
        & (x1-x3).*(y-y3)-(x-x3).*(y1-y3)>=-tol*sqrt((x1-x3).^2+(y1-y3).^2),1);
    if ~isempty(FindIndex)
        TriangleIndex(i)=index(FindIndex);
    end
    
%    FindIndex=index(FindIndex);
%    TriangleIndex(i)=FindIndex(1);
%     TriangleIndex(i)=index(find(FindIndex,1));
    
    
end


FindIndex=find(~isnan(TriangleIndex));
if ~isempty(FindIndex)
    num=length(FindIndex);
    three_p=3*num;
    index1=1:3:three_p;
    index2=2:3:three_p;
    index3=3:3:three_p;
    t1=T(TriangleIndex(FindIndex),1);
    t2=T(TriangleIndex(FindIndex),2);
    t3=T(TriangleIndex(FindIndex),3);
    M=sparse([reshape(repmat(index1,3,1),1,[]),repmat(index2,1,3),repmat(index3,1,3)],...
        [1:three_p,repmat([index1,index2,index3],1,2)],...
        [ones(three_p,1);V(t1,1);V(t2,1);V(t3,1);V(t1,2);V(t2,2);V(t3,2)]);
    Bary(FindIndex,:)=reshape(M\reshape([ones(num,1),X(FindIndex),Y(FindIndex)]',[],1),3,[])';
    
end


% %tic
% [row,column]=size(X);
% if row>column
%     X=X';
% end
% 
% [row,column]=size(Y);
% if row>column
%     Y=Y';
% end
% 
% n=length(X);
% m=size(T,1);
% TriangleIndex=zeros(n,1);
% 
% for i=1:n
% X1=repmat(X(i)-tol,m,1);
% X2=repmat(X(i)+tol,m,1);
% Y1=repmat(Y(i)-tol,m,1);
% Y2=repmat(Y(i)+tol,m,1);
% index=find(Analyze(:,1)<=X2 & Analyze(:,2)>=X1 & Analyze(:,3)<=Y2 & Analyze(:,4)>=Y1);
% 
% 
% x=X(i);
% y=Y(i);
% 
% 
% x1=V(T(mod(index-1,m)+1,1),1);
% y1=V(T(mod(index-1,m)+1,1),2);
% x2=V(T(mod(index-1,m)+1,2),1);
% y2=V(T(mod(index-1,m)+1,2),2);
% x3=V(T(mod(index-1,m)+1,3),1);
% y3=V(T(mod(index-1,m)+1,3),2);
% 
% % the three sides of the triangle
% FindIndex=((x2-x1).*(y3-y1)-(y2-y1).*(x3-x1)).*((x2-x1).*(y-y1)-(y2-y1).*(x-x1))>=-tol & ((x3-x2).*(y1-y2)-(y3-y2).*(x1-x2)).*((x3-x2).*(y-y2)-(y3-y2).*(x-x2))>=-tol & ((x3-x1).*(y2-y1)-(y3-y1).*(x2-x1)).*((x3-x1).*(y-y1)-(y3-y1).*(x-x1))>=-tol;
% 
% 
% TriangleIndex(i)=mod(index(FindIndex)-1,m)+1;
% 
% end
% 
% %toc
% 
% 

return;

% Below is old code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[row,column]=size(X);
if row>column
    X=X';
end

[row,column]=size(Y);
if row>column
    Y=Y';
end

n=length(X);
m=size(T,1);

X1=reshape(repmat(X-tol,m,1),[],1);
X2=reshape(repmat(X+tol,m,1),[],1);
Y1=reshape(repmat(Y-tol,m,1),[],1);
Y2=reshape(repmat(Y+tol,m,1),[],1);
Analyze=repmat(Analyze,n,1);
index=find(Analyze(:,1)<=X2 & Analyze(:,2)>=X1 & Analyze(:,3)<=Y2 & Analyze(:,4)>=Y1);


x=X(floor((index-1)/m)+1)';
y=Y(floor((index-1)/m)+1)';

[row,column]=size(x);
if row<column
    x=x';
    y=y';
end

x1=V(T(mod(index-1,m)+1,1),1);
y1=V(T(mod(index-1,m)+1,1),2);
x2=V(T(mod(index-1,m)+1,2),1);
y2=V(T(mod(index-1,m)+1,2),2);
x3=V(T(mod(index-1,m)+1,3),1);
y3=V(T(mod(index-1,m)+1,3),2);

% the three sides of the triangle
FindIndex=((x2-x1).*(y3-y1)-(y2-y1).*(x3-x1)).*((x2-x1).*(y-y1)-(y2-y1).*(x-x1))>=-tol & ((x3-x2).*(y1-y2)-(y3-y2).*(x1-x2)).*((x3-x2).*(y-y2)-(y3-y2).*(x-x2))>=-tol & ((x3-x1).*(y2-y1)-(y3-y1).*(x2-x1)).*((x3-x1).*(y-y1)-(y3-y1).*(x-x1))>=-tol;

TriangleIndex=zeros(n,1);
TriangleIndex(floor((index(FindIndex)-1)/m)+1)=mod(index(FindIndex)-1,m)+1;

% for i=1:n
%     x=X(i);
%     y=Y(i);
%     x1=x+tol;
%     x2=x-tol;
%     y1=y+tol;
%     y2=y-tol;
%     index1=find(Analyze(:,1)<=x1 & Analyze(:,2)>=x2);
%     index2=find(Analyze(:,3)<=y1 & Analyze(:,4)>=y2);
%     index3=ismember(index1,index2);
%     index4=find(index3);
% 
%     
%     m=size(index4);
%     for j=1:m
%         tri_index=index1(index4(j));
%         x1=V(T(tri_index,1),1);
%         y1=V(T(tri_index,1),2);
%         x2=V(T(tri_index,2),1);
%         y2=V(T(tri_index,2),2);
%         x3=V(T(tri_index,3),1);
%         y3=V(T(tri_index,3),2);
%         
%         % the first side of the triangle
%         if (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)>0
%             if (x2-x1)*(y-y1)-(y2-y1)*(x-x1)+tol<=0
%                 continue;
%             end
%         elseif (x2-x1)*(y-y1)-(y2-y1)*(x-x1)-tol>=0
%             continue;
%         end
%         
%         % the second side of the triangle
%         if (x3-x2)*(y1-y2)-(y3-y2)*(x1-x2)>0
%             if (x3-x2)*(y-y2)-(y3-y2)*(x-x2)+tol<=0
%                 continue;
%             end
%         elseif (x3-x2)*(y-y2)-(y3-y2)*(x-x2)-tol>=0
%             continue;
%         end
%         
%         % the third side of the triangle
%         if (x3-x1)*(y2-y1)-(y3-y1)*(x2-x1)>0
%             if (x3-x1)*(y-y1)-(y3-y1)*(x-x1)+tol<=0
%                 continue;
%             end
%         elseif (x3-x1)*(y-y1)-(y3-y1)*(x-x1)-tol>=0
%             continue;
%         end
%         
%         TriangleIndex(i)=tri_index;
%         break;
%         
%     end
% end




end

