%%% ECON 634 Macro II
%%% Problem Set 2
%%% solve the NGM with tech shock using VFI
%%% Yaobin Wang
%%% 9/19/2017

%%% Question 1
%%% Functional Equation (FE):
%%% v(A_i,k) = max{u[A_i*k^alpha+(1-delta)k-k']+beta*E[v(A',k')]}
%%% where u(c) = c^(1-sigma)/(1-sigma), return function
%%%       E[v(A',k'}] = pi_ih*v_h(A_h',k')+pi_il*v_l(A_l',k')
%%%       i = h, l
%%% State variables: k, A
%%% Control variable: k', c

%%% Question 2
clear all;
close all;
clc;

%%%% Set up parameters
alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;

pi_hh = 0.977;
pi_ll = 0.926;
PI = [pi_hh 1-pi_hh; 1-pi_ll pi_ll]; %transition matrix

A_h = 1.1;
A_l = 0.678;

%%%% Set up discretized state space
k_min = 0;
k_max = 45;
num_k = 1000; % number of points in the grid for k

k = linspace(k_min, k_max, num_k);

k_mat = repmat(k', [1 num_k]); % this will be useful in a bit

%%%% Set up consumption and return function
% consumption function
% 1st dim(rows): k today, 2nd dim (cols): k' chosen for tomorrow
cons_h = A_h*(k_mat .^ alpha) + (1 - delta) * k_mat - k_mat'; 
cons_l = A_l*(k_mat .^ alpha) + (1 - delta) * k_mat - k_mat'; 

% return function
ret_h = cons_h .^ (1 - sigma) / (1 - sigma); 
ret_l = cons_l .^ (1 - sigma) / (1 - sigma); 
% negative consumption is not possible -> make it irrelevant by assigning
% it very large negative utility
ret_h(cons_h < 0) = -Inf;
ret_l(cons_l < 0) = -Inf;

%%%% Iteration
dis_h = 1;
dis_l = 1;
tol = 1e-06; % tolerance for stopping 
v_h_guess = zeros(1, num_k);
v_l_guess = zeros(1, num_k);
n = 1; % loop counter
tic; % record the time elapsed
while dis_h > tol || dis_l > tol
    % compute the utility value for all possible combinations of A, k amd k'
    value_mat_h = ret_h + beta * (pi_hh*repmat(v_h_guess, [num_k 1]) + ...
        (1-pi_hh)*repmat(v_l_guess, [num_k 1]));
    value_mat_l = ret_l + beta * ((1-pi_ll)*repmat(v_h_guess, [num_k 1])...
        + pi_ll*repmat(v_l_guess, [num_k 1]));
    
    % find the optimal k' for every k and A;
    [vfn_h, pol_indx_h] = max(value_mat_h, [], 2);
    vfn_h = vfn_h';
    [vfn_l, pol_indx_l] = max(value_mat_l, [], 2);
    vfn_l = vfn_l';
    
    % what is the distance between current guess and value function
    dis_h = max(abs(vfn_h - v_h_guess));
    dis_l = max(abs(vfn_l - v_l_guess));
    
    % if distance is larger than tolerance, update current guess and
    % continue, otherwise exit the loop
    v_h_guess = vfn_h;
    v_l_guess = vfn_l;
    
    % count the number of loops
    n = n+1;
end
time = toc;

% report time and loop counter
fprintf('Time to solve the value function iteration was %2.2f seconds.\n',time)
fprintf('It took %4.0f loops to solve the value function iteration. \n',n)

% policy function
g_h = k(pol_indx_h); 
g_l = k(pol_indx_l); 

% savings
sav_h = g_h - (1-delta)*k;
sav_l = g_l - (1-delta)*k;

% plot the value function over K for each state of A
figure
plot(k,vfn_h,k,vfn_l)
legend('state A^h','state A^l','Location','SouthEast')
xlabel('k')
ylabel('value function')
title('The value function over K for A^h and A^l')
% The value function is increasing and concave over K for each state of A.

%%% Question 3
% plot the policy function over K for each state of A
figure
plot(k,g_h,k,g_l)
legend('state A^h','state A^l','Location','SouthEast')
xlabel('k')
ylabel('policy function (k'')')
title('The policy function (K'') over K for A^h and A^l')
% The policy function is increasing in K and A.

% plot savings over K for each state of A
figure
plot(k,sav_h,k,sav_l)
legend('state A^h','state A^l','Location','NorthWest')
xlabel('k')
ylabel('savings')
title('Savings over K for A^h and A^l')
% Savings is increasing in A but not monotonic in K. Savings for each state
% of A are firstly increasing and then decreasing in K.

%%% Question 4
%%% Simulation with A_h = 1.1 and A_l = 0.678
t = 1000;
X = rand(t-1,1);
A = zeros(t,1);
A(1) = A_h; %assign A_0
i = 1;
while i < t
    if A(i) == A_h
        if X(i) < pi_hh
            A(i+1) = A_h;
        else
            A(i+1) = A_l;
        end
    else
        if X(i) < pi_ll
            A(i+1) = A_l;
        else
            A(i+1) = A_h;
        end
    end
    i = i+1;
end

K = zeros(t+1,1);
k_index = 2;
K(1) = k(k_index); %assign k_0
i = 1;
while i < t
    if A(i) == A_h
        K(i+1) = g_h(k_index);
    else
        K(i+1) = g_l(k_index);
    end
    k_index = find(k == K(i+1));
    i = i+1;
end

Y = zeros(t,1);
i = 1;
while i < t+1;
    Y(i) = A(i)*(K(i)^alpha);
    i = i+1;
end

std_Y = std(Y(100:t));