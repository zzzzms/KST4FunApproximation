%This is a final version for constructing KB spline called LKB to
%approximate 2D functions in a list of testing functions. It is written by
%Dr. Ming-Jun Lai based on the 2D KST construction from Zhaiming Shen who
%is a Ph.D. student under Dr. Ming-Jun Lai's supervision in 2018--2023. 
%This code is divided into two parts.  
% The first part is to generate the
%least squares matrix X_data based on 101^2 equally-spaced points over [0, 1] and
% generate LKB splines once for all.  
%Then we use the X_data and find a discrete least squares fit to each testing function 
%and compute the accuracy of LKB spline approximation of various testing functions.
%The second part is to find the best data locations based on a matrix cross 
% approximation method developed by Dr. Ming-Jun Lai with a help from Mr. Kenneth Allen 
% who was a Ph.D. studet in 2018--2021 and solve a much small
% size of discrete least squares fit. The RMSEs are computed and compared with the 
% standard discrete least squares fit in the first part.  

%you can comment off phi_q after the first run
phi_q;
% number n of linear spline basis and n = length(LKB)
n = 101; h = 2/(n-1); xx = linspace(0,2,n);
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

% compute LKB basis functions by denoising  KB functions to get LKB functions.  
% and then generate a least squares matrix.
for i =1:n
[c,d,V,T,Analyze,tol] = LaisSplines4denoisingPos(xx,yy,KB{i}(xx,yy));
nVal = SplineEvaluation2D(V, T, Analyze, c, d, xx(:), yy(:), tol);
LKB{i}=reshape(nVal,size(xx)); 
LKB2{i}=c; %The  spline coefficient vector
end

%The least squares matrix
X_data = zeros((1/hh+1)^2, n);
for i = 1:n
    X_data(:,i) = reshape(LKB{i},[],1);
end

%The following is for evaluation. 
hhh = 0.0025;
xxx = linspace(0,1,1/hhh+1);yyy = linspace(0,1,1/hhh+1);
[xxx,yyy] = meshgrid(xxx,yyy);

%Compute a standard discrete least square fit. 
KBSE=zeros(20,4);nn=zeros(20,1); %error matrix.
for k=1:20
ny=testfunctions2d(xx,yy,k); 
y_LS = reshape(ny(xx,yy),[],1);
coeff = lsqminnorm(X_data, y_LS); %coeff to be sent via internet.
I=find(abs(coeff)>=1e-8);nn(k)=size(I,1);
f_reconst = 0; %after received the coeff from the internet. 
for i =1:n
    f_reconst = f_reconst + coeff(i)*LKB2{i};
end
% check errors at the given data locations
spVal = SplineEvaluation2D(V,T,Analyze,f_reconst, d, xxx(:), yyy(:), tol);
exact=ny(xxx,yyy); 
%figure, subplot(131), SplV=reshape(spVal,401,401); surf(xxx,yyy,SplV), shading interp
%subplot(132), surf(xxx,yyy,exact), shading interp
%subplot(133), surf(xxx,yyy,exact-SplV), shading interp
RMSE=rms((exact(:)-spVal)./(max(1,abs(exact(:)))));
Rel_L2_error = norm(exact(:)-spVal)/norm(exact(:));
Linf=norm(exact(:)-spVal,inf)/norm(exact(:),inf);
KBSE(k,:)=[k,RMSE,Rel_L2_error,Linf];
end

%The following is part 2;
Phi=X_data; Phi=full(Phi);  err=1e-4; 
I=isnan(Phi); J=find(I==1); Phi(J)=0; 
[m,L]=size(Phi); A=zeros(L,1); NNr=zeros(20,1);
for i=1:L
    a(i)=norm(Phi(:,i)); 
end
Iind=find(abs(a)>=1e-6); 
nPhi=Phi(:,Iind); [n1,r]=size(nPhi);

%This is the part we use a matrix cross approximation.
Im=LaisCrossApproximation(nPhi);

nPhi2=nPhi(Im,:); 
 NNr=ones(20,1); NP=size(Im,1);
KBSME=zeros(20,4);
for k=1:20
ny=testfunctions2d(xx,yy,k); 
y_LS = reshape(ny(xx(Im),yy(Im)),[],1);
%x0= myOMP(nPhi2,y_LS,500); 
%x0=nPhi2\y_LS;
x0 = lsqminnorm(nPhi2, y_LS);
ny_LS=reshape(ny(xx,yy),[],1);
norm(nPhi*x0-ny_LS,inf)
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
b=KBSME(:,2)./KBSE(:,2); 
figure, subplot(122), plot(b,'r')
title(['The ratio of the accuracies based on ',num2str(NP),' vs 10201 sampled values'])
hold on
plot(b,'ro')
xlabel('testing function labels')
subplot(121), plot(KBSE(:,2),'bo')
title('RMSEs of testing functions')
return
