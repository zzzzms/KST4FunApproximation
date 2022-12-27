function c= SplineInterpolation2D( V, T, Analyze, d, r, X, Y, Z, lambda, tol, caseN, opts )
% This matlab code is written by Yidong Xu under Dr. Ming-Jun Lai's
% supervision in 2019.  
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
switch caseN
   case 1
        c=ConstrainedLeastSquareKKT(A,Z,K,H,lambda,tol);
   case 2
%         c=IterativeProjection(A,Z,K,H,lambda,NIterations,tol);
        c=ConstrainedLeastSquareGradientProjection(A,Z,K,H,lambda,NIterations,tol);
    case 3
        c=SlowDecompositionDantzigWolfeSolveEquations(V,T,d,A,Z,K,H,lambda,x_partition_num,y_partition_num,Boundary,NIterations,Anticycle,InfNum,tol);
    case 4
        c=AnotherDecompositionDantzigWolfeSolveEquations(V,T,d,A,Z,K,H,lambda,x_partition_num,y_partition_num,Boundary,NIterations,Anticycle,InfNum,tol);
    case 5
        c=DecompositionProjection(V,T,d,A,Z,K,H,lambda,NIterations,x_partition_num,y_partition_num,tol);
end
warning('on','MATLAB:rankDeficientMatrix');
end


function [ x ] = ConstrainedLeastSquareKKT( A, b, K, H, lambda, tol )
% [row1,column1]=size(A);
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
% M=[H',A'*A+lambda*K;sparse(row2,row2),H;sparse(1,column2+row2)];
% N=[A'*sparse(b);sparse(row2+1,1)];
% x=M\N;
% x=x(row2+1:row2+column1);

% "t" means transpose.
% n=size(A,2);
% AtA=A'*A+spdiags(1e-12*ones(n,1),0,n,n);
% AtA=A'*A;
% Atb=A'*b;
% AtA_inv_Atb=AtA\Atb;
% AtA_inv_Ht=AtA\H';
% 
% lambda=(H*AtA_inv_Ht)\(H*AtA_inv_Atb);
% 
% x=AtA\(Atb-H'*lambda);
% x(1:6,:)
end






function [ x ] = ConstrainedLeastSquareGradientProjection( A, b, K, H, lambda, NIterations, tol )

b=sparse(b);
[row,column]=size(H);
if isempty(H)
    x=[A'*A;sparse(1,column)]\[A'*b;0];
    return;
end

AtA_plus_lambdaK=A'*A+lambda*K;
Atb=A'*b;

t=1/max(eig(AtA_plus_lambdaK));
H0=sparse(row,1);
x=sparse(column,1);
HHT=H*H';

for i=1:NIterations
    x_new=ProjectionOntoAffineSpaceDirectly(x-t*(AtA_plus_lambdaK*x-Atb),H,H0,HHT,tol);
    if norm(x-x_new,inf)<tol
        x=x_new;
        break;
    end
    x=x_new;
end

disp('Number of iterations=');
disp(i);

if i==NIterations
    disp('Reached the max number of iterations.');
end

end







function [ x ] = ProjectionOntoAffineSpaceDirectly( x0, A, b, AAT, tol )
% min||x-x0||^2 such that Ax=b

[row,~]=size(A);
x0=sparse(x0);

% x=x0-A'*([A*A';sparse(1,row)]\[A*x0-sparse(b);0]);
x=x0-A'*([AAT;sparse(1,row)]\[A*x0-sparse(b);0]);

if norm(A*x-b,inf)>tol
    disp('The solution may not be correct.');
end


end



function [ x ] = SlowDecompositionDantzigWolfeSolveEquations( V, T, d, A, b, K, H, lambda, ...
    x_partition_num, y_partition_num, Boundary, NIterations, Anticycle, InfNum, tol )

[L,DA,N,column_permutation]=Decomposition2D(V,T,d,A,b,K,H,lambda,x_partition_num,y_partition_num,tol);
x=SlowDantzigWolfeSolveEquations(L,DA,N,Boundary,NIterations,Anticycle,InfNum,tol);
[~,inverse]=sort(column_permutation);
x=x(inverse);
x=x(1:size(A,2));

end






function [ x ] = AnotherDecompositionDantzigWolfeSolveEquations( V, T, d, A, b, K, H, lambda, ...
    x_partition_num, y_partition_num, Boundary, NIterations, Anticycle, InfNum, tol )

[L,DA,N,column_permutation]=Decomposition2D(V,T,d,A,b,K,H,lambda,x_partition_num,y_partition_num,tol);
x=AnotherDantzigWolfeSolveEquations(L,DA,N,Boundary,NIterations,Anticycle,InfNum,tol);
[~,inverse]=sort(column_permutation);
x=x(inverse);
x=x(1:size(A,2));

end









function [ x ] = DecompositionProjection( V, T, d, A, b, K, H, lambda, NIterations, x_partition_num, y_partition_num, tol )

[L,DA,N,column_permutation]=...
    Decomposition2D(V,T,d,A,b,K,H,lambda,x_partition_num,y_partition_num,tol);
x=AlternatingProjection(L,DA,N,NIterations,tol);
[~,inverse]=sort(column_permutation);
x=x(inverse);
x=x(1:size(A,2));
end



function [ x ] =AlternatingProjection(L,A,b,NIterations,tol)
b=sparse(b); m=size(L{1,1},1); N=size(A,1);
BlockSize=zeros(N,2);
BlockRowColumnIndices=cell(N,2);
index1=m; index2=0;
L_combined=[];
for i=1:N
    [row,column]=size(A{i,1});
    BlockSize(i,:)=size(A{i,1});
    BlockRowColumnIndices{i,1}=index1+1:index1+row;
    index1=index1+row;
    BlockRowColumnIndices{i,2}=index2+1:index2+column;
    index2=index2+column;
    L_combined=[L_combined,L{i,1}];
end

n=sum(BlockSize(:,2));
first_m_indices=1:m;
x=sparse(n,1);
LLT=L_combined*L_combined';
AAT=cell(N,1);
for i=1:N
    AAT{i,1}=A{i,1}*A{i,1}';
end

for i=1:NIterations
    x_new=ProjectionOntoAffineSpaceDirectly(x,L_combined,b(first_m_indices),LLT,tol);
    for j=1:N
    x_new(BlockRowColumnIndices{j,2})=ProjectionOntoAffineSpaceDirectly...
       (x_new(BlockRowColumnIndices{j,2}),A{j,1},b(BlockRowColumnIndices{j,1}),AAT{j,1},tol);
    end
    if norm(x-x_new,inf)<tol
        x=x_new;
        break;
    end
    x=x_new;
end

disp('Number of iterations=');
disp(i);
if i==NIterations
    disp('Reached the max number of iterations.');
end
end




function [ x ] = AlternatingProjection_Test( L, A, b, NIterations, tol )

b=sparse(b);
m=size(L{1,1},1);
N=size(A,1);

BlockSize=zeros(N,2);
BlockRowColumnIndices=cell(N,2);
index1=m;
index2=0;
L_combined=[];
for i=1:N
    [row,column]=size(A{i,1});
    BlockSize(i,:)=size(A{i,1});
    BlockRowColumnIndices{i,1}=index1+1:index1+row;
    index1=index1+row;
    BlockRowColumnIndices{i,2}=index2+1:index2+column;
    index2=index2+column;
    L_combined=[L_combined,L{i,1}];
end

n=sum(BlockSize(:,2));
first_m_indices=1:m;

x=sparse(n,1);
x_new=x;
y=x;
y_new=x;
d_new=zeros(N,1);
d=d_new;
LLT=L_combined*L_combined';
AAT=cell(N,1);
for i=1:N
    AAT{i,1}=A{i,1}*A{i,1}';
end

for i=1:2
    x_new=ProjectionOntoAffineSpaceDirectly(y,L_combined,b(first_m_indices),LLT,tol);
    for j=1:N
        y_new(BlockRowColumnIndices{j,2})=ProjectionOntoAffineSpaceDirectly...
            (x_new(BlockRowColumnIndices{j,2}),A{j,1},b(BlockRowColumnIndices{j,1}),AAT{j,1},tol);
        d_new(j)=norm(x_new(BlockRowColumnIndices{j,2})-y_new(BlockRowColumnIndices{j,2}));
    end
    x=x_new;
    y=y_new;
    d=d_new;
end

for i=1:NIterations
%     perm=randperm(N+1);
%     for j=1:N+1
%         k=perm(j);
%         if k>N
%             x_new=ProjectionOntoAffineSpaceDirectly(x_new,L_combined,b(first_m_indices),LLT,tol);
%         else
%             x_new(BlockRowColumnIndices{k,2})=ProjectionOntoAffineSpaceDirectly...
%                 (x_new(BlockRowColumnIndices{k,2}),A{k,1},b(BlockRowColumnIndices{k,1}),AAT{k,1},tol);
%         end
%         
%     end
    
%     x_new=ProjectionOntoAffineSpaceDirectly(x,L_combined,b(first_m_indices),LLT,tol);
%     for j=1:N
%         x_new(BlockRowColumnIndices{j,2})=ProjectionOntoAffineSpaceDirectly...
%             (x_new(BlockRowColumnIndices{j,2}),A{j,1},b(BlockRowColumnIndices{j,1}),AAT{j,1},tol);
%     end
    
    x_new=ProjectionOntoAffineSpaceDirectly(y,L_combined,b(first_m_indices),LLT,tol);
    for j=1:N
        y_new(BlockRowColumnIndices{j,2})=ProjectionOntoAffineSpaceDirectly...
            (x_new(BlockRowColumnIndices{j,2}),A{j,1},b(BlockRowColumnIndices{j,1}),AAT{j,1},tol);
        d_new(j)=norm(x_new(BlockRowColumnIndices{j,2})-y_new(BlockRowColumnIndices{j,2}));
        y_new(BlockRowColumnIndices{j,2})=y(BlockRowColumnIndices{j,2})+...
            (d(j)/(d(j)-d_new(j)))*(y_new(BlockRowColumnIndices{j,2})-y(BlockRowColumnIndices{j,2}));
        
    end
    
    
    
%     [maax,indexx]=max(abs(x-x_new))
%     if norm(x-x_new,inf)<tol
%         x=x_new;
%         break;
%     end
    norm(y-y_new,inf)
    if norm(y-y_new,inf)<tol
        x=y_new;
        break;
    end
    x=x_new;
    y=y_new;
    d=d_new;
    
%     x=x_new;
end

disp('Number of iterations=');
disp(i);

if i==NIterations
    disp('Reached the max number of iterations.');
end

end



% function [ x ] = IterativeProjection( A, b, K, H, lambda, NIterations, tol )
% 
% [row,column]=size(H);
% if isempty(H)
%     x=[A'*A;sparse(1,column)]\[A'*sparse(b);0];
%     return;
% end
% 
% AtA_plus_lambdaK=A'*A+lambda*K;
% Atb=A'*b;
% % AtA=A;
% % Atb=b;
% H0=sparse(row,1);
% x0=sparse(column,1);
% x=x0;
% 
% for i=1:NIterations
%     x0=ProjectionOntoAffineSpace(x0,AtA_plus_lambdaK,Atb,tol);
%     x1=ProjectionOntoAffineSpace(x0,H,H0,tol);
%     
% %     if norm(x0-x1,inf)<tol
%     if norm(x-x1,inf)<tol
%         x=x1;
%         break;
%     end
%     
%     x0=x1;
%     x=x1;
%     
% end
% 
% disp('Number of iterations=');
% disp(i);
% 
% if i==NIterations
%     disp('Reached the max number of iterations.');
% end
% 
% 
% 
% end






% function [ x ] = ProjectionOntoAffineSpace( x0, A, b, tol )
% % min||x-x0||^2 such that Ax=b
% 
% [row,column]=size(A);
% M=[speye(column),A';A,sparse(row,row);sparse(1,column+row)];
% N=[x0;b;sparse(1,1)];
% x=M\N;
% 
% if norm(M*x-N,inf)>tol
%     disp('The solution may not be correct.');
% end
% 
% x=x(1:column);
% 
% 
% 
% end

