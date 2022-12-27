Lambda = 1/sqrt(2);

% choose k_max = 2 or 3 or 4
%k_max = 2;
%k_max = 3;
k_max = 4;

[alpha_0, beta_0, mapping_0] = phi_construct(k_max,0);
[alpha_1, beta_1, mapping_1] = phi_construct(k_max,1);
[alpha_2, beta_2, mapping_2] = phi_construct(k_max,2);
[alpha_3, beta_3, mapping_3] = phi_construct(k_max,3);
[alpha_4, beta_4, mapping_4] = phi_construct(k_max,4);

Alpha = {alpha_0, alpha_1, alpha_2, alpha_3, alpha_4};
Beta = {beta_0, beta_1, beta_2, beta_3, beta_4};
Mapping = {mapping_0, mapping_1, mapping_2, mapping_3, mapping_4};

x_0 = [];
y_0 = [];
x_1 = [];
y_1 = [];
x_2 = [];
y_2 = [];
x_3 = [];
y_3 = [];
x_4 = [];
y_4 = [];


for point_x = keys(mapping_0)
    x_0 = [x_0, point_x{1}];
    y_0 = [y_0, mapping_0(point_x{1})];
end
for point_x = keys(mapping_1)
    x_1 = [x_1, point_x{1}];
    y_1 = [y_1, mapping_1(point_x{1})];
end
for point_x = keys(mapping_2)
    x_2 = [x_2, point_x{1}];
    y_2 = [y_2, mapping_2(point_x{1})];
end
for point_x = keys(mapping_3)
    x_3 = [x_3, point_x{1}];
    y_3 = [y_3, mapping_3(point_x{1})];
end
for point_x = keys(mapping_4)
    x_4 = [x_4, point_x{1}];
    y_4 = [y_4, mapping_4(point_x{1})];
end


max([min(x_0), min(x_1), min(x_2), min(x_3), min(x_4)])
min([max(x_0), max(x_1), max(x_2), max(x_3), max(x_4)])


% phi_q functions
phi_0 = @(x) interp1(x_0,y_0,x);
%phi_0 = @(x) spline(x_0,y_0,x);
phi_1 = @(x) interp1(x_1,y_1,x);
%phi_1 = @(x) spline(x_1,y_1,x);
phi_2 = @(x) interp1(x_2,y_2,x);
%phi_2 = @(x) spline(x_2,y_2,x);
phi_3 = @(x) interp1(x_3,y_3,x);
%phi_3 = @(x) spline(x_3,y_3,x);
phi_4 = @(x) interp1(x_4,y_4,x);
%phi_4 = @(x) spline(x_4,y_4,x);







