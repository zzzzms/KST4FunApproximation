function [LKB,LKB2]=LKBsplines(KB)
%This function generates  denoised KB splines
%from the input KB splines. It is written by 
% Dr. Ming-Jun Lai on May, 2022. 
% number of data locations over [0, 1]^2. %hh = 1/(length(LKB{1})-1); %%
hh = 0.01;
xx = linspace(0,1,1/hh+1);yy = linspace(0,1,1/hh+1);
[xx,yy] = meshgrid(xx,yy);
n=size(KB,1);
% compute LKB basis functions by denoising  KB functions to get LKB functions.  
% and then generate a least squares matrix.
for i =1:n
[c,d,V,T,Analyze,tol] = LaisSplines4denoisingPos(xx,yy,KB{i}(xx,yy));
nVal = SplineEvaluation2D(V,T,Analyze,c,d,xx(:), yy(:), tol);
LKB{i}=reshape(nVal,size(xx)); 
LKB2{i}=c; %The  spline coefficient vector
end