function [ V, T ] = Triangulation2D( V, B, H, C, d, tol )
%[ V, T ] = Triangulation2D( V, B, H, C, d, tol )
% Written by Yidong Xu under Dr. Ming-Jun Lai's supervision in Fall, 2018.
% V=[x1,y1; x2,y2; ...] vertices
% B=[1;2;3;...] vertex indices of the boundary (counterclockwise, no repetition)
% H={[1;2;...];[3;4;...];...} vertex indices of holes (clockwise, no repetition)
% C=[1,2; 3,4; ...] vertex indices of constrained edges

m1=length(B);
m2=0;
for i=1:size(H,1)
    m2=m2+length(H{i,1});
end
m3=size(C,1);
m=m1+m2+m3;
E=zeros(m,2);
E(1:m1,:)=[B,[B(2:m1);B(1)]];
% for i=1:m1-1
%     E(i,:)=[B(i),B(i+1)];
% end
% E(m1,:)=[B(m1),B(1)];
index=m1;
for i=1:size(H,1)
    s=length(H{i,1});
    E(index+1:index+s,:)=[H{i,1},[H{i,1}(2:s);H{i,1}(1)]];
    index=index+s;
end
E(m1+m2+1:m,:)=C;



polygon=[V(B,:);V(B(1),:);zeros(m2+2*size(H,1),2)];
index=m1+1;
for i=1:size(H,1)
    s=length(H{i,1});
    polygon(index+1:index+s+2,:)=[[NaN,NaN];V(H{i,1},:);V(H{i,1}(1),:)];
    index=index+s+2;
end
% plot(polygon(:,1),polygon(:,2),'LineWidth',2);



% Pick vertices in the grid.
x_min=min(V(:,1));
x_max=max(V(:,1));
y_min=min(V(:,2));
y_max=max(V(:,2));

new_x=x_min:d:x_max;
new_y=y_min:d:y_max;
[new_x,new_y]=meshgrid(new_x,new_y);
% new_x=reshape(new_x,[],1);
% new_y=reshape(new_y,[],1);
new_x=new_x(:);
new_y=new_y(:);

% figure,plot(new_x,new_y,'r*');
% Delete those vertices outside the polygon (possibly nonconvex).
in=InPolygon2D(new_x,new_y,polygon(:,1),polygon(:,2),tol);
new_x(~in)=[];
new_y(~in)=[];

% Delete those vertices too close to the edges. ('too close' means
% distance<=d/2)
distance=d/2;
for i=1:m
    v1=E(i,1);
    v2=E(i,2);
    x1=V(v1,1);
    y1=V(v1,2);
    x2=V(v2,1);
    y2=V(v2,2);
    FindIndex1=find(new_x>=min(x1,x2)-distance-tol & new_x<=max(x1,x2)+distance+tol...
        & new_y>=min(y1,y2)-distance-tol & new_y<=max(y1,y2)+distance+tol);
    X=new_x(FindIndex1);
    Y=new_y(FindIndex1);
    
    distance2=zeros(length(FindIndex1),1);
    FindIndex2=sign((y2-y1)*(Y-y1)+(x2-x1)*(X-x1)).*sign((y2-y1)*(Y-y2)+(x2-x1)*(X-x2))<=0;
    distance2(FindIndex2)=abs((x2-x1)*(Y(FindIndex2)-y1)-(X(FindIndex2)-x1)*(y2-y1))/sqrt((x2-x1)^2+(y2-y1)^2);
    distance2(~FindIndex2)=min([sqrt((X(~FindIndex2)-x1).^2+(Y(~FindIndex2)-y1).^2),...
        sqrt((X(~FindIndex2)-x2).^2+(Y(~FindIndex2)-y2).^2)],[],2);
    FindIndex2=distance2<=distance+tol;
    
    FindIndex2=FindIndex1(FindIndex2);
    
    new_x(FindIndex2)=[];
    new_y(FindIndex2)=[];
end
V=[V;new_x,new_y];


% Pick uniformly spaced vertices on each edge.
% Preallocating will keep getting larger.
preallocating_vertices=1000;
num_of_V=size(V,1);
V=[V;zeros(preallocating_vertices,2)];
preallocating_splittededges=preallocating_vertices;
num_of_splittededges=0;
splittedE=zeros(preallocating_splittededges,2);

for i=1:m
    v1=E(i,1);
    v2=E(i,2);
    x1=V(v1,1);
    y1=V(v1,2);
    x2=V(v2,1);
    y2=V(v2,2);
    distance=sqrt((x2-x1)^2+(y2-y1)^2);
    num=floor(distance/d+tol)-1;
    if num>0
        coef=(1:num)';
        points=((num+1-coef)*[x1,y1]+coef*[x2,y2])/(num+1);
        
        if num_of_splittededges+num+1>size(splittedE,1)
            splittedE=[splittedE;zeros(num+1+preallocating_splittededges,2)];
            preallocating_splittededges=preallocating_splittededges*2;
        end
        splittedE(num_of_splittededges+1:num_of_splittededges+num+1,:)...
            =[[v1;(num_of_V+1:num_of_V+num)'],[(num_of_V+1:num_of_V+num)';v2]];
        num_of_splittededges=num_of_splittededges+num+1;
        
        if num_of_V+num>size(V,1)
            V=[V;zeros(num+preallocating_vertices,2)];
            preallocating_vertices=preallocating_vertices*2;
        end
        V(num_of_V+1:num_of_V+num,:)=points;
        num_of_V=num_of_V+num;
        
    else
        if num_of_splittededges+1>size(splittedE,1)
            splittedE=[splittedE;zeros(1+preallocating_splittededges,2)];
            preallocating_splittededges=preallocating_splittededges*2;
        end
        splittedE(num_of_splittededges+1,:)=[v1,v2];
        num_of_splittededges=num_of_splittededges+1;
        
    end
    
end


V(num_of_V+1:end,:)=[];
% V=V(1:num_of_V,:);

DT=delaunayTriangulation(V(:,1),V(:,2),splittedE(1:num_of_splittededges,:));
V=DT.Points;
T=DT.ConnectivityList;

% Now delete unqualified triangles.
% incenters=incenter(DT,(1:size(DT,1))');
incenters=incenter(DT);
in=InPolygon2D(incenters(:,1),incenters(:,2),polygon(:,1),polygon(:,2),tol);
T(~in,:)=[];


% Delete degenerate triangles.
X=[V(T(:,1),1),V(T(:,2),1),V(T(:,3),1)];
Y=[V(T(:,1),2),V(T(:,2),2),V(T(:,3),2)];
areas=abs(((X(:,2)-X(:,1)).*(Y(:,3)-Y(:,1))-(X(:,3)-X(:,1)).*(Y(:,2)-Y(:,1)))/2);
FindIndex1=areas<=tol;
T(FindIndex1,:)=[];

end

