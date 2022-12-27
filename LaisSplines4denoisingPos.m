function [nc, d, V, T, Analyze, tol] = LaisSplines4denoisingPos(xx,yy,zz)
%This is to demo how to use bivariate splines to denoise a noised surface. 
%It is written by Dr. Ming-Jun Lai on May 28, 2019. It is modified by Dr.
%Ming-Jun Lai to preserve the nonnegativity.  
%close all
warning('off');
%you can adjust the following parameters.
r=2; d=3*r+2; tol=1e-5; lambda=1;
%zz=myK2(xx,yy)/5;
X=xx(:); Y=yy(:);  Z=zz(:);
%figure, surf(xx,yy,zz), shading interp, title('A Noised Surface')
%return
opts=struct('NIterations',0,'x_partition_num',0,'y_partition_num',0,'Boundary',0,'Anticycle',0,'InfNum',0);
opts.NIterations=1000;
caseN=1; %number for solution methods.
opts.x_partition_num=3;opts.y_partition_num=3;
P=[0 0;1 0;1 1;0 1]; B=[1;2;3;4]; H=[]; C=[]; dist=0.1; tol=1e-6;
[V, T] = Triangulation2D(P, B, H, C, dist, tol );
Initialization2D(d);  tol=1e-2;
[T,Analyze]= AnalyzeTriangulation2D(V,T);
c = SplineInterpolation2D(V,T,Analyze,d,r,X,Y,Z,lambda,tol,caseN,opts);
I=find(c>=-1); nc=zeros(size(c)); nc(I)=c(I);
