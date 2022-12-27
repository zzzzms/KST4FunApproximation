function [I,k] = simple_greedy_maxvol(X,I_initial)
%finds a close to dominant rxr submatrix of nxr matrix X
%swaps up to two of the largest elements each itteration
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
A = X(I,:); %sub-matrix

%if abs(det(A)) < 1e-10
%    error('initial submatrix is not inverible')
%end

for k=1:1000
    Y = X/A;
    Y_abs = abs(Y);
    
    [y,linear_index] = max(Y_abs(:)); %y = largest in magnitude entry in Y
    if y < 1+epsilon
        break
    elseif k == 1000
        disp('simple_greedy_maxvol did not converge in 1000 steps.')
    end
    [i1,j1] = ind2sub([n,r],linear_index); %converts linear index to matrix index
    I(j1) = i1; %replace j1th row with i1th row
    
    Y_abs(:,j1) = 0; %makes sure next largest value is not in the same column
    Y_abs(i1,:) = 0; %makes sure next largest value is not in the same row
    [~,linear_index] = max(Y_abs(:)); %finds largest entry not in column i1 or row j1
    [i2,j2] = ind2sub([n,r],linear_index); %converts linear index to matrix index
    
    if abs(det(Y([i1 i2],[j1 j2]))) > y %makes sure volume is increasing
        I(j2) = i2; %replaces jth row with ith row
    end

    A = X(I,:);
end