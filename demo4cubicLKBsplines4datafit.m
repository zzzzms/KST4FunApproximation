%This is a demo code for using cubic LKB splines to approximate smooth global functions. 
%It is based on LKB-splines which are obtained by replacing the K-outer function 
% by cubic B-splines to have KB splines and then by denoising the KB splines.
% K-outer functions are from the Kolmogorov representation theorem.  
% This demo code is written by Dr. Ming-Jun Lai based on the 2D KST 
% construction from Zhaiming Shen who is a Ph.D. student under Dr. Ming-Jun Lai's 
% supervision in 2018--2023. 
% The first part of the demo is to generate the KB and LKB splines.
% Then the code is divided into two parts for two data fitting strategies.  
% The first strategy is to solve a discrete least squares fit of various functions based
% on 101x101 equally-spaced points over [0, 1]^2 and the function values.  
% The second strategy is to generate a set of magic locations with the set containing 
% much less than 101x101 points. Based on the magic data locations, we solve the 
% discrete least squares fit for the same testing functions in the first strategy. 
% Finally we compare the accuracies of LKB spline approximation of various testing 
% functions based on 101x101 points vs based the magic set. 
% The comparison shows that the accuracies are similar, the ratios of the
% accuracies from two strategies are between 1 to 3.

% The following line is the beginning of this demo code. 
phi_q;%you can comment off phi_q after the first run

% number n of linear spline basis and n = length(LKB)
n = 201; M=100; 

%We first generate cubic KB splines
KB=CubicKBsplines(n,Lambda,phi_0,phi_1,phi_2,phi_3,phi_4);

% number of data locations over [0, 1]^2. %hh = 1/(length(LKB{1})-1); %%
hh = 0.01; xx = linspace(0,1,1/hh+1);yy = linspace(0,1,1/hh+1);
[xx,yy] = meshgrid(xx,yy);
LaisSplines4denoising %this generates cubic LKB splines and save them in LKB and LKB2. 

%The discrete least squares matrix based on 101x101 equally-spaced
%points over [0, 1]^2.
X_data = zeros((1/hh+1)^2, n);
for i = 1:n
    X_data(:,i) = reshape(LKB{i},[],1);
end

hhh = 0.0025; %We will use this set of points for evaluation.
xxx = linspace(0,1,1/hhh+1);yyy = linspace(0,1,1/hhh+1);
[xxx,yyy] = meshgrid(xxx,yyy);

%strategy 1, we use a compressive sensing approach to find the sparse solution. 
KBCS=zeros(M,5); %error matrix. 
Phi=X_data; Phi=full(Phi);  %err=1e-4; 
I=isnan(Phi); J=find(I==1); Phi(J)=0; 
[m,L]=size(Phi); A=zeros(L,1); NN=zeros(30,1);
for i=1:L
    a=norm(Phi(:,i)); Phi(:,i)=Phi(:,i)/a; %normalization
    A(i)=a;
end
I=isnan(Phi); J=find(I==1); Phi(J)=0; 
for k=1:M
ny=testfunctions_2d(xx,yy,k);
y_LS = reshape(ny(xx,yy),[],1);
x0= myOMP(Phi,y_LS,500); norm(Phi*x0-y_LS,inf);
I=find(abs(x0)>=1e-8); [L, size(I)], NN(k)=size(I,1);
coeff=zeros(L,1); coeff(I)=x0(I)./A(I); %coeff to be sent via internet. 
f_reconst = 0;
for i =1:n
    f_reconst = f_reconst + coeff(i)*LKB2{i};
end
spVal = SplineEvaluation2D(V,T,Analyze,f_reconst, d, xxx(:), yyy(:), tol);
exact=ny(xxx,yyy); 
%figure, subplot(131), SplV=reshape(spVal,401,401); surf(xxx,yyy,SplV), shading interp
%subplot(132), surf(xxx,yyy,exact), shading interp
%subplot(133), surf(xxx,yyy,exact-SplV), shading interp
RMSE=rms((exact(:)-spVal(:))./max(1,abs(exact(:))));
Rel_L2_error = norm(exact(:)-spVal(:))/norm(exact(:));
Linf=norm(exact(:)-spVal(:),inf)/norm(exact(:),inf);
KBCS(k,:)=[k,RMSE,Rel_L2_error,Linf NN(k)];
end

%strategy 2: We first find a magic set.
Phi=X_data; Phi=full(Phi);  err=1e-4; 
I=isnan(Phi); J=find(I==1); Phi(J)=0; 
[m,L]=size(Phi); A=zeros(L,1); NNr=zeros(20,1);
for i=1:L
    a(i)=norm(Phi(:,i)); 
end
Iind=find(abs(a)>=1e-6); 
nPhi=Phi(:,Iind); [n1,r]=size(nPhi);

%The following computation is based on a matrix cross approximation.
Im=LaisCrossApproximation(nPhi);
figure, plot(xx(Im),yy(Im),'r*')
title('The magic locations for data fittinig based on KST')
nPhi2=nPhi(Im,:); 
 NNr=ones(M,1); NP=size(Im',1);
KBSME=zeros(M,4);
for k=1:M
ny=testfunctions_2d(xx,yy,k); 
y_LS = reshape(ny(xx(Im),yy(Im)),[],1);
%x0 = lsqminnorm(nPhi2, y_LS);
x0=myOMP(nPhi2,y_LS,500);
coeff=zeros(L,1); coeff(Iind)=x0;
I=find(abs(coeff)>=1e-8);NNr(k)=size(I,1);
f_reconst = 0;
for i =1:n
 f_reconst = f_reconst + coeff(i)*LKB2{i};
end
spVal = SplineEvaluation2D(V,T,Analyze,f_reconst, d, xxx(:), yyy(:), tol);
exact=ny(xxx,yyy); 
%figure, subplot(131), SplV=reshape(spVal,401,401); surf(xxx,yyy,SplV), shading interp
%subplot(132), surf(xxx,yyy,exact), shading interp
%subplot(133), surf(xxx,yyy,exact-SplV), shading interp
RMSE=rms((exact(:)-spVal)./max(1,abs(exact(:))));
Rel_L2_error = norm(exact(:)-spVal)/norm(exact(:));
Linf=norm(exact(:)-spVal,inf)/norm(exact(:),inf);
KBSME(k,:)=[k,RMSE,Rel_L2_error,Linf];
end
%[KBCS(M,2), KBSME(M,2)]
%return
% M=106;
% b1=KBSME(1:M,2);c1=KBCS(1:M,2); 
% b2=KBSME(1:M,3);c2=KBCS(1:M,3);  
% b3=KBSME(1:M,4);c3=KBCS(1:M,4);
% figure, subplot(131), hold on, 
% plot(c1,'r'), plot(b1,'b'),
% legend('RMSE101^2','RMSE99')
% subplot(132),plot(c2,'r'), hold on
% plot(b2,'b'),legend('reL2-101^2','rel2-99')
% subplot(133), plot(c3,'r'), hold on
% plot(b3,'b'),  legend('rel-Linf101^2','rel-Linf99')
% axis tight

% figure, subplot(121), plot(xx(Im),yy(Im),'r*'), title('A Magic Data Locations')
% subplot(122), hold on,plot(c1,'r'), plot(b1,'b'),
% legend('RMSE101^2','RMSE99'),
% title(['RMSEs of DLS from 101^2 vs ',num2str(NP), ' points'])
M=100;
b1=KBSME(1:M,2);c1=KBCS(1:M,2); 
figure, subplot(121), plot(xx(Im),yy(Im),'r*'), grid on
title('A Set of Magic Data Locations')
subplot(122), hold on,plot(c1,'r'), plot(b1,'b'),
legend('RMSE101^2','RMSEmagic'),
title(['RMSEs of DLS from 101^2 vs ',num2str(NP), ' points'])
return

