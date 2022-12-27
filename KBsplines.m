function KB=KBsplines(n)
%This generates KBsplines based on linear B-splines.
%It is written by Zhaiming Shen under Dr. Ming-Jun Lai
%supervision in 2021.
%e.g. n = 101; 
h = 2/(n-1); xx = linspace(0,2,n);
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

