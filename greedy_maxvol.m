function [I,k] = greedy_maxvol(X,I_initial)
%finds a close to dominant rxr submatrix of nxr matrix X
%swaps up to r of the largest elements each itteration
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

I = I_initial; %positions of initial sub-matrix
A = X(I,:); %initial sub-matrix. Can change to any initial condition

%if abs(det(A)) < 1e-10
%    error('initial submatrix is not inverible')
%end

i = zeros(r,1); %row indices
j = zeros(r,1); %column indices

for k=1:1000
    Y = X/A;
    Y_abs = abs(Y);
    
    [y,linear_index] = max(Y(:)); %y = largest entry
    if y <= 1+epsilon
        break
    elseif k == 1000
        disp('greedy_maxvol did not converge in 1000 steps.')
    end
    [i(1),j(1)] = ind2sub([n,r],linear_index); %index of y
    I(j(1)) = i(1); %replace jth row with ith row
    Y_abs(:,j(1)) = 0; %makes sure next largest value is not in the same column
    Y_abs(i(1),:) = 0; %makes sure next largest value is not in the same row
    
    for s=2:r
        [~,linear_index] = max(Y_abs(:)); %searches for largest index not in the original column
        [i(s),j(s)] = ind2sub([n,r],linear_index); %[i,j] = index
        
        V = Y(i(1:(s-1)),j(1:(s-1)));
        B = Y(i(1:(s-1)),j(s));
        C = Y(i(s),j(1:(s-1)));
        y = Y(i(s),j(s));
        if  abs(y-C*(V\B)) > 1 %makes sure volume is increasing
            I(j(s)) = i(s); %replaces jth row with ith row
            Y_abs(:,j(s)) = 0;
            Y_abs(i(s),:) = 0;
        else
            break
        end
    end
    A = X(I,:);
end