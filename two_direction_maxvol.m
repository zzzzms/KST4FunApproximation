function [I,J,num_iter] = two_direction_maxvol(X,I_initial,J_initial)
%two directional maxvol
%gives indices of an rxr submatrix of mxn with maximum determinant in modulus
%A = X(I,J)
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

num_iter = 0;

for k=1:1000
    Y = X(:,J)/A;
    Z = A\X(I,:);
    num_iter = num_iter+2;
    [y,linear_index_I] = max(abs(Y(:))); %Ilinear = linear index
    [z,linear_index_J] = max(abs(Z(:))); %Ilinear = linear index
    if max([y z]) < 1+epsilon
        break
    elseif k==1000
        error('did not converge in 1000 steps')
    end
    
    if y > z
        [i,j] = ind2sub([m,r],linear_index_I); %[i,j] = index of largest entry in columns
        I(j) = i; %replace jth row with ith row
    else
        [p,q] = ind2sub([r,n],linear_index_J); %[p,q] = index of largest entry in rows
        J(p) = q;
    end
    A = X(I,J);
end