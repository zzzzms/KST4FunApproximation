function [I,k] = maxvol(X,I_initial)
%finds a close to dominant rxr submatrix of nxr matrix X
%swaps the largest element each itteration
%A = X(I,:)
%k = # of steps to converge
%It is written by Kenneth Allen under Dr. Ming-Jun Lai's supervision in
%2020.

epsilon = 1e-8;
[n,r] = size(X);
[r_check1,r_check2] = size(I_initial);
if max(r_check1,r_check2) ~= r
    error('The size of the initial indices must equal the width of the matrix')
end

I = I_initial; %index set of sub-matrix
A = X(I,:); %initial sub-matrix

%if abs(det(A)) < 1e-12
%    error('initial submatrix is not inverible')
%end

for k=1:1000
    Y = abs(X/A);
    [y,linear_index] = max(Y(:));
    if y <= 1+epsilon
        break
    elseif k==1000
        disp('maxvol did not converge in 1000 steps.')
    end
    
    [i,j] = ind2sub([n,r],linear_index); %converts linear index to matrix index
    I(j) = i; %replace jth row with ith row
    A = X(I,:);
end