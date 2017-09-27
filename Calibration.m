%%% ECON 634 Macro II
%%% Problem Set 2
%%% solve the NGM with tech shock using VFI
%%% Yaobin Wang
%%% 9/19/2017

%%% Question 4
%%% Calibrating A_h and A_l inside the model
clear all;
close all;
clc;

alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;
    
pi_hh = 0.977;
pi_ll = 0.926;
PI = [pi_hh 1-pi_hh; 1-pi_ll pi_ll];
PI_con = PI^100;

k_min = 0;
k_max = 45;
num_k = 1000;
k = linspace(k_min, k_max, num_k);
k_mat = repmat(k', [1 num_k]);
    
A_h = 1.002; %start with a value close to 1 to boost the calculation time
A_l = (1- (PI_con(1,1)*A_h))/PI_con(1,2); %to ensure that the averge of A is 1

dis_std = 1;
tol_std = 3e-03;
tic;

while dis_std > tol_std
    cons_h = A_h*(k_mat .^ alpha) + (1 - delta) * k_mat - k_mat'; 
    cons_l = A_l*(k_mat .^ alpha) + (1 - delta) * k_mat - k_mat'; 

    ret_h = cons_h .^ (1 - sigma) / (1 - sigma); 
    ret_l = cons_l .^ (1 - sigma) / (1 - sigma); 
    ret_h(cons_h < 0) = -Inf;
    ret_l(cons_l < 0) = -Inf;
    
    dis_h = 1;
    dis_l = 1;
    tol = 1e-06;
    v_h_guess = zeros(1, num_k);
    v_l_guess = zeros(1, num_k);
    
    while dis_h > tol || dis_l > tol
        value_mat_h = ret_h + beta * (pi_hh*repmat(v_h_guess, [num_k 1]) + ...
            (1-pi_hh)*repmat(v_l_guess, [num_k 1]));
        value_mat_l = ret_l + beta * ((1-pi_ll)*repmat(v_h_guess, [num_k 1])...
            + pi_ll*repmat(v_l_guess, [num_k 1]));
    

        [vfn_h, pol_indx_h] = max(value_mat_h, [], 2);
        vfn_h = vfn_h';
        [vfn_l, pol_indx_l] = max(value_mat_l, [], 2);
        vfn_l = vfn_l';
    
        dis_h = max(abs(vfn_h - v_h_guess));
        dis_l = max(abs(vfn_l - v_l_guess));
    
        v_h_guess = vfn_h;
        v_l_guess = vfn_l;
    end

    g_h = k(pol_indx_h); 
    g_l = k(pol_indx_l); 

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

    std_Y = std(Y(100:t)); %cut off first 100 observations

    dis_std = std_Y - 0.018;
    
    A_h = A_h - 1e-04;
    A_l = (1- PI_con(1,1)*A_h)/PI_con(1,2);
end
time = toc;
fprintf('Time to calibrate A_h and A_l was %3.2f seconds.\n',time)

A_h = A_h + 1e-04;
A_l = (1- PI_con(1,1)*A_h)/PI_con(1,2);

display (A_h);
display (A_l);
display (std_Y);

% Multiple calibration results imply that for A_h between 1 to 1.0013 and
% for A_l between 0.9958 to 1, respectively, the standard deviation of
% output (Y) is around 2%.



