%This is to demo how to use bivariate splines to denoise a noised surface. 
%It is written by Dr. Ming-Jun Lai on May 28, 2019. It is modified by Dr.
%Ming-Jun Lai to denoise the surface of KB splines.  
warning('off');
r=2; d=3*r+2; tol=1e-5; lambda=1; %you can adjust the following parameters.
X=xx(:); Y=yy(:);  
P=[0 0;1 0;1 1;0 1]; B=[1;2;3;4]; H=[]; C=[]; dist=0.1; tol=1e-6;
[V, T] = Triangulation2D(P, B, H, C, dist, tol );
Initialization2D(d); 
[T,Analyze]= AnalyzeTriangulation2D(V,T);
SplineInterpolation2Dpart1
n=size(KB,2); X=xx(:);Y=yy(:);
LKB = cell(1,n); LKB2 = cell(1,n);
for k=1:n
  W=KB{k}(xx,yy); W=W(:);
%c = SplineInterpolation2D(V,T,Analyze,d,r,X,Y,Z,lambda,tol,caseN,opts);
%nVal = SplineEvaluation2D(V, T, Analyze, c, d, xx(:), yy(:), tol);
%I=find(c3>=-1); nc=zeros(size(c3)); nc(I)=c3(I);
c3=(dA)\(A'*sparse(W)); 
LKB2{k}=c3/200000;
nVal3 = SplineEvaluation2D(V,T,Analyze,c3,d,xx(:),yy(:),tol);
LKB{k}=reshape(nVal3,size(xx))/200000; 
end