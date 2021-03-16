function [peak_max, out] = peak_lin_sys(Acl, F, X0, order, signed, rank_tol, Tmax)
%PEAK_LIN_SYS Find the maximum absolute value of any element of a query 
%F*x for a linear system x' = Acl x, where x starts in the convex hull of
%initial points X0. (closed-loop if controller is applied, so Acl = A - B K)
%Holds for arbitrary switches x' = Acl_i x for subsystems {Acl_i}
%
% Input:
%   Acl:        Linear System x' = Acl x (or cell for switching)
%   F:          Query matrix, f_j = Fj'x along trajectories
%   X0:         Set of initial feasible points
%   order:      Order of moment relaxation (degree=2*order)
%   signed:     max (Fjx) (signed=1) or |Fjx| (signed=0, default)
%   ranktol:    rank tolerance of moment matrix 
%   Tmax:       Maximum time. By default is Infinite (time-independent).
%               Tmax ~= Inf will recover optimal time at optimality, and
%               allow for time-varying safety contours
%Ouputs:
%   peak_max:   Maximum absolute value of query abs(Fj x) along
%               trajectories starting in X0
%   out:        Data structure that contains the peak estimation information


%% Parse the input
if nargin < 7
    Tmax = Inf;    
end

if nargin < 6
    rank_tol = 1e-3;
end

if nargin < 5
    signed = 0;
end

if nargin < 4
    order = 2;
end
       
d = order*2;


if ~iscell(Acl)
    Acl = {Acl};
end

nx0 = size(X0, 2);
nf = size(F, 1);

nsys = length(Acl);

%% Solve the problem

opt = cell(nx0, nf);
peak_val = zeros(nx0, nf);


for i = 1:nx0
    X0i = X0(:, i);          
    for j = 1:nf
        
        Fj = F(j, :);  
        %form an equivalent impulse-response problem
        
        Bcl = X0i;
        Ccl = Fj;
        
        %run peak code on I/O pair (i, j)
        [peak_val_curr, opt_curr] = peak_impulse_siso(Acl, Bcl, Ccl, order, signed, rank_tol, Tmax);
        peak_val(i, j) = peak_val_curr;
        opt{i, j} = opt_curr;
    end
end

[peak_max, ind_max] = max(peak_val, [], 1);

opt_max = cell(nf, 1);
optimal = zeros(nf, 1);
time_max = zeros(nf, 1);
for j = 1:nf
    
    opt_max{j} = opt{ind_max(j), j};
    optimal(j) = opt_max{j}.optimal;
    time_max(j) = opt_max{j}.tp;
end

%package to output
out = struct;
out.peak_val = peak_val;
out.peak_max = peak_max;
out.opt = opt;
out.opt_max = opt_max;
out.optimal = optimal;
out.time_max = time_max;




end

