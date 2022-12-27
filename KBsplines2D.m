%This is a demo code for using KST to approximate smooth global functions. 
%It is based on KB-splines. That is, replacing the K-outer function by B-splines,
% the Kolmogorov representation produces a continuous function which is called 
% KB-spline.  As these KB-splines are not smooth, bivariate splines are used 
% to smooth these KB splines to produce LKB splines. 
% These LKB splines are used to approximate 2D functions in a list of testing 
% functions. This demo code is written by Dr. Ming-Jun Lai based on the 2D KST 
% construction from Zhaiming Shen who is a Ph.D. student under Dr. Ming-Jun Lai's 
% supervision in 2018--2023. 
% This code is divided into two parts.  The first part is to generate the
%LKB splines and a least squares matrix once for all.  The second part is to
%test the accuracy of LKB spline approximation of various testing functions.
%Instead of the standard least squares approach, we adopt the sparse
%solution technique (myOMP.m) to have a minimum number of LKB splines for
%each testing function.  
%you can comment off phi_q after the first run
phi_q;
% number n of linear spline basis and n = length(LKB)
n = 21; h = 2/(n-1); xx = linspace(0,2,n);
% construct linear splines B
B = cell(1,n);KB = cell(1,n); LKB = cell(1,n); LKB2 = cell(1,n); %%
for i =1:n
    if i==1
        B{i} = @(x) (1/h)*(h-x).*(x<=h & x>=0);
        KB{i} = @(x,y) B{i}(phi_0(x)+Lambda*phi_0(y))+B{i}(phi_1(x)+Lambda*phi_1(y))+B{i}(phi_2(x)+Lambda*phi_2(y))+B{i}(phi_3(x)+Lambda*phi_3(y))+B{i}(phi_4(x)+Lambda*phi_4(y));
    elseif i == n
        B{i} = @(x) (1/h)*(x-2+h).*(x<=2 & x>=(2-h));
        KB{i} = @(x,y) B{i}(phi_0(x)+Lambda*phi_0(y))+B{i}(phi_1(x)+Lambda*phi_1(y))+B{i}(phi_2(x)+Lambda*phi_2(y))+B{i}(phi_3(x)+Lambda*phi_3(y))+B{i}(phi_4(x)+Lambda*phi_4(y));
    else
        B{i} = @(x) (1/h)*(x-xx(i)+h).*(x<=xx(i) & x>xx(i)-h) + (1/h)*(xx(i)+h-x).*(xx(i)<x & x<=xx(i)+h);
        KB{i} = @(x,y) B{i}(phi_0(x)+Lambda*phi_0(y))+B{i}(phi_1(x)+Lambda*phi_1(y))+B{i}(phi_2(x)+Lambda*phi_2(y))+B{i}(phi_3(x)+Lambda*phi_3(y))+B{i}(phi_4(x)+Lambda*phi_4(y));
    end
end

% number of data locations over [0, 1]^2. %hh = 1/(length(LKB{1})-1); %%
hh = 0.01;
xx = linspace(0,1,1/hh+1);yy = linspace(0,1,1/hh+1);
[xx,yy] = meshgrid(xx,yy);

I=zeros(size(xx));
for i=1:n
    Val=KB{i}(xx,yy); I=I+Val; J=find(Val<-1e-1); size(J), end
mesh(xx,yy,I)

% compute LKB basis functions and generate a least squares matrix.
% for i =1:n
% [c, d, V, T, Analyze, tol] = LaisSplines4denoisingPos(xx,yy,KB{i}(xx,yy));
% nVal = SplineEvaluation2D(V, T, Analyze, c, d, xx(:), yy(:), tol);
% LKB{i}=reshape(nVal,size(xx)); 
% LKB2{i}=c; %The  spline coefficient vector
% end
LaisSplines4denoising %generate LKB splines and save them in LKB and LKB2. 

%The least squares matrix
X_data = zeros((1/hh+1)^2, n);
for i = 1:n
    X_data(:,i) = reshape(LKB{i},[],1);
end
%the above is part 1. the following is  part2
hhh = 0.0025;
xxx = linspace(0,1,1/hhh+1);yyy = linspace(0,1,1/hhh+1);
[xxx,yyy] = meshgrid(xxx,yyy);
%Compute a least square fit. 
% KBSE=zeros(10,4); %error matrix.
% for k=1:20
% ny=testfunctions_2d(xx,yy,k); 
% y_LS = reshape(ny(xx,yy),[],1);
% coeff = lsqminnorm(X_data, y_LS); %coeff to be sent via internet.
% size(coeff)
% f_reconst = 0; %after received the coeff from the internet. 
% for i =1:n
%     f_reconst = f_reconst + coeff(i)*LKB2{i};
% end
% % check errors at the given data locations
% spVal = SplineEvaluation2D(V,T,Analyze,f_reconst, d, xxx(:), yyy(:), tol);
% exact=ny(xxx,yyy); 
% %figure, subplot(131), SplV=reshape(spVal,401,401); surf(xxx,yyy,SplV), shading interp
% %subplot(132), surf(xxx,yyy,exact), shading interp
% %subplot(133), surf(xxx,yyy,exact-SplV), shading interp
% RMSE=rms(exact(:)-spVal);
% Rel_L2_error = norm(exact(:)-spVal)/norm(exact(:));
% Linf=norm(exact(:)-spVal,inf)/norm(exact(:),inf);
% KBSE(k,:)=[k,RMSE,Rel_L2_error,Linf];
% end

%Nex use the compressive sensing approach to find the sparse solution. 
KBCS=zeros(20,4); %error matrix. 
Phi=X_data; Phi=full(Phi);  err=1e-4; 
I=isnan(Phi); J=find(I==1); Phi(J)=0; 
[m,L]=size(Phi); A=zeros(L,1); NN=zeros(10,1);
for i=1:L
    a=norm(Phi(:,i)); Phi(:,i)=Phi(:,i)/a; %normalization
    A(i)=a;
end
I=isnan(Phi); J=find(I==1); Phi(J)=0; 
for k=1:20
ny=testfunctions_2d(xx,yy,k); 
y_LS = reshape(ny(xx,yy),[],1);
%coeff = lsqminnorm(X_data, y_LS); %coeff to be sent via internet.
x0= myOMP(Phi,y_LS,500); norm(Phi*x0-y_LS,inf)
I=find(abs(x0)>=1e-8); [L, size(I)], NN(k)=size(I,1); 
coeff=zeros(L,1); coeff(I)=x0(I)./A(I); 
f_reconst = 0;
for i =1:n
    f_reconst = f_reconst + coeff(i)*LKB2{i};
end
spVal = SplineEvaluation2D(V,T,Analyze,f_reconst, d, xxx(:), yyy(:), tol);
exact=ny(xxx,yyy); 
%figure, subplot(131), SplV=reshape(spVal,401,401); surf(xxx,yyy,SplV), shading interp
%subplot(132), surf(xxx,yyy,exact), shading interp
%subplot(133), surf(xxx,yyy,exact-SplV), shading interp
RMSE=rms(exact(:)-spVal);  %./max(1,abs(exact(:))));
Rel_L2_error = norm(exact(:)-spVal)/norm(exact(:));
%Linf=norm(exact(:)-spVal,inf); %);
Linf=norm(exact(:)-spVal,inf)/norm(exact(:),inf);
KBCS(k,:)=[k,RMSE,Rel_L2_error,Linf];
end
%Compare the solution from LSQ and Sparse Solutions.
KBCS


return
for i=1:n
A=LKB{i}; figure, subplot(1,2,1), surf(xx,yy,A), shading interp
b=KB{i}; B=b(xx,yy); subplot(1,2,2), surf(xx,yy,B), shading interp
end

return
B=imread('bank.pgm'); B=double(B);
ny=B(100:1:300,100:1:300); 
y_LS=ny(:); x1=Phi\y_LS; norm(Phi*x1- y_LS,inf)
x0= OMP(Phi,y_LS,500); norm(Phi*x0-y_LS,inf)

return