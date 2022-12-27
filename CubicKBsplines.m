function KB=CubicKBsplines(n,Lambda,phi_0,phi_1,phi_2,phi_3,phi_4)
%This generates KB splines based on cubic B-spline functions.
%n is the number of cubic spline basis. 
%This is based on the code from Zhaiming Shen under the direction of Dr.
%Ming-Jun Lai in 2022. 
%n = 21; %for example.  
h = 2/(n-1); xx = linspace(0,2,n);
% construct C^2 cubic B splines 
B = cell(1,n+2);KB = cell(1,n+2); LKB = cell(1,n+2); LKB2 = cell(1,n+2); %%
B{1} = @(x) CubicSplines(x,xx(1)-3*h,xx(2)-3*h,xx(3)-3*h,xx(4)-3*h,xx(5)-3*h);
B{2} = @(x) CubicSplines(x,xx(1)-2*h,xx(2)-2*h,xx(3)-2*h,xx(4)-2*h,xx(5)-2*h);
B{3} = @(x) CubicSplines(x,xx(1)-1*h,xx(2)-1*h,xx(3)-1*h,xx(4)-1*h,xx(5)-1*h);
B{n} = @(x) CubicSplines(x,xx(n-4)+h,xx(n-3)+h,xx(n-2)+h,xx(n-1)+h,xx(n)+h);
B{n+1} = @(x) CubicSplines(x,xx(n-4)+2*h,xx(n-3)+2*h,xx(n-2)+2*h,xx(n-1)+2*h,xx(n)+2*h);
B{n+2} = @(x) CubicSplines(x,xx(n-4)+3*h,xx(n-3)+3*h,xx(n-2)+3*h,xx(n-1)+3*h,xx(n)+3*h);

KB{1} = @(x,y) B{1}(phi_0(x)+Lambda*phi_0(y))+B{1}(phi_1(x)+Lambda*phi_1(y))+B{1}(phi_2(x)+Lambda*phi_2(y))+B{1}(phi_3(x)+Lambda*phi_3(y))+B{1}(phi_4(x)+Lambda*phi_4(y));   
KB{2} = @(x,y) B{2}(phi_0(x)+Lambda*phi_0(y))+B{2}(phi_1(x)+Lambda*phi_1(y))+B{2}(phi_2(x)+Lambda*phi_2(y))+B{2}(phi_3(x)+Lambda*phi_3(y))+B{2}(phi_4(x)+Lambda*phi_4(y));   
KB{3} = @(x,y) B{3}(phi_0(x)+Lambda*phi_0(y))+B{3}(phi_1(x)+Lambda*phi_1(y))+B{3}(phi_2(x)+Lambda*phi_2(y))+B{3}(phi_3(x)+Lambda*phi_3(y))+B{3}(phi_4(x)+Lambda*phi_4(y));   
KB{n} = @(x,y) B{n}(phi_0(x)+Lambda*phi_0(y))+B{n}(phi_1(x)+Lambda*phi_1(y))+B{n}(phi_2(x)+Lambda*phi_2(y))+B{n}(phi_3(x)+Lambda*phi_3(y))+B{n}(phi_4(x)+Lambda*phi_4(y));   
KB{n+1} = @(x,y) B{n+1}(phi_0(x)+Lambda*phi_0(y))+B{n+1}(phi_1(x)+Lambda*phi_1(y))+B{n+1}(phi_2(x)+Lambda*phi_2(y))+B{n+1}(phi_3(x)+Lambda*phi_3(y))+B{n+1}(phi_4(x)+Lambda*phi_4(y));   
KB{n+2} = @(x,y) B{n+2}(phi_0(x)+Lambda*phi_0(y))+B{n+2}(phi_1(x)+Lambda*phi_1(y))+B{n+2}(phi_2(x)+Lambda*phi_2(y))+B{n+2}(phi_3(x)+Lambda*phi_3(y))+B{n+2}(phi_4(x)+Lambda*phi_4(y));   

for i=4:n-1
    B{i} = @(x) CubicSplines(x,xx(i-3),xx(i-2),xx(i-1),xx(i),xx(i+1));
    KB{i} = @(x,y) B{i}(phi_0(x)+Lambda*phi_0(y))+B{i}(phi_1(x)+Lambda*phi_1(y))+B{i}(phi_2(x)+Lambda*phi_2(y))+B{i}(phi_3(x)+Lambda*phi_3(y))+B{i}(phi_4(x)+Lambda*phi_4(y));   
end


