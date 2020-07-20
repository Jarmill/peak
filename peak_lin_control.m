function [peak_val, out] = peak_lin_control(A, B, K, X0, order, rank_tol, Tmax)
%PEAK_LIN_CONTROL Find the maximum control effort for a state-feedback
%controller u = -Kx for arbitrary switching of the linear systems:
% x' = Ai x + Bi u = (Ai - Bi K) x. Assume that initial set of points is 
%polytopic with a corners in X0.
%
% Input:
%   A,B:        Linear System
%   order:      Order of moment relaxation (degree=2*order)
%   ranktol:    rank tolerance of moment matrixof 
%   Tmax:       Maximum time. By default is Infinite (time-independent).
%               Tmax ~= Inf will recover optimal time at optimality, and
%               allow for time-varying safety contours
%Ouputs:
%   peak_val:   Maximum absolute value of control effort
%   out:        Data structure that contains the peak control information



if nargin < 6
    Tmax = Inf;    
end

if nargin < 5
    rank_tol = 1e-3;
end

if nargin < 4
    order = 2;
end

if ~iscell(A)
    A = {A};
    B = {B};
end

%nx = size(A{1}, 1);
nu = size(B{1}, 2);
nx0 = size(X0, 2);

nsys = length(A);
Acl = cell(nsys, 1); %closed loop dynamics
        
d = order*2;

opt = cell(nx0, nu);
peak_val = zeros(nx0, nu);

for i = 1:nx0
    X0i = X0(:, i);          
    for j = 1:nu
        
        Kj = K(j, :);  
        %form an equivalent impulse-response problem
        for k = 1:nsys
            Acl{k} = A{k} - B{k}*K;
        end
        
        Bcl = X0i;
        Ccl = Kj;
        
        %run peak code on I/O pair (i, j)
        [peak_val_curr, opt_curr] = peak_impulse_siso(Acl, Bcl, Ccl, order);
        peak_val(i, j) = peak_val_curr;
        opt{i, j} = opt_curr;
    end
end

[peak_val, ind_max] = max(peak_val, [], 1);

opt_max = cell(nu, 1);
optimal = zeros(nu, 1);
for j = 1:nu
    
    opt_max{j} = opt{ind_max(j), j};
    optimal(j) = opt_max{j}.optimal;
end

%package to output
out = struct;
out.peak_val = peak_val;
out.opt = opt;
out.opt_max = opt_max;
out.optimal = optimal;





end

