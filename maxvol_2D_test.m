%test for the maxvol and alt_maxvol algorithms
%D1 = importdata('solps_trial_data.mat'); % 1D data
%M1 = D1.one_d_low_noise'; %low noise
%[n,r] = size(M1);

m = 200; %mxn
n = 200;
r = 50; %rank
N = 1; %number of random matrices to generate

t_2_av = 0; %maxvol average time
t_alt_av = 0; %simple greedy maxvol average time
t_galt_av = 0; %greedy maxvol average time

k_2_av = 0; %maxvol averate  # iterations
k_alt_av = 0; %simple greedy maxvol average # iterations
k_galt_av = 0; %greedy maxvol average # iterations

for i=1:N
X = rand(m,n);
I_initial = randperm(m,r);
J_initial = randperm(n,r);

tic
[I_2,J_2,k_2] = two_direction_maxvol(X,I_initial,J_initial);
t_2 = toc; %time to compute maxvol

t_2_av = t_2_av+t_2;
k_2_av = k_2_av+k_2;

tic
[I_alt,J_alt,k_alt] = alt_maxvol(X,I_initial,J_initial);
t_alt = toc; %time to compute simple greedy maxvol

t_alt_av = t_alt_av+t_alt;
k_alt_av = k_alt_av+k_alt;

tic
[I_galt,J_galt,k_galt] = alt_greedy_maxvol(X,I_initial,J_initial);
tg = toc; %time to compute greedy maxvol

t_galt_av = t_galt_av+tg;
k_galt_av = k_galt_av+k_galt;
end

%log10(abs(det(X(I_initial,J_initial))))
%log10(abs(det(X(I_2,J_2))))
%log10(abs(det(X(I_alt,J_alt))))
%log10(abs(det(X(I_galt,J_galt))))

t_2_av/N
t_alt_av/N
t_galt_av/N

k_2_av/N
k_alt_av/N
k_galt_av/N