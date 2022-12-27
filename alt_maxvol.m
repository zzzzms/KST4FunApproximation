function [I,J,num_iter] = alt_maxvol(X,I_initial,J_initial)
%alternating maxvol
%gives indices of an rxr submatrix of mxn with maximum determinant in modulus
%num_iter = number of backslash operations to converge

epsilon = 1e-8; %error tolerance
[m,n] = size(X);

I = I_initial; %index set of sub-matrix rows
J = J_initial; %index set of sub-matrix columns
A = X(I,J); %initial sub-matrix

r = size(I);
r = max(r);
r2 = size(J);
r2 = max(r2);
if r~=r2
    error('initial submatrix is not rxr.')
end

row_dom = 0; %indicates if near dominant in rows
num_iter = 0; %number of iterations to converge

if cond(A) > 1e12
    error('Initial A is not invertible')
end

for k=1:1000
    Y = abs(X(:,J)/A);
    num_iter = num_iter+1;
    [y,linear_index_I] = max(Y(:)); %y is max in columns
    if y > 1+epsilon %if the volume change by swapping one column
        [i,j] = ind2sub([m r],linear_index_I); %[i,j] = index
        I(j) = i; %replace jth row with ith row
        A = X(I,J);
        column_dom = 0; %indicates not near dominant in columns
    elseif row_dom == 1
        break %if near dominant in rows and columns
    else
        column_dom = 1; %indicates near dominant in columns
    end
    
    Z = abs(A\X(I,:));
    num_iter = num_iter+1;
    [z,linear_index_J] = max(Z(:)); %z is max in rows
    if z > 1+epsilon
        [p,q] = ind2sub([r n],linear_index_J); %[i,j] = index
        J(p) = q; %replace ith column with jth row
        A = X(I,J);
        row_dom = 0; %indicates not near dominant in rows
    elseif column_dom == 1
        break %if near dominant in rows and columns
    else
        row_dom = 1; %indicates near dominant in rows
    end
    
    if k==1000
        error('alt_maxvol did not converge in 1000 steps')
    end
end