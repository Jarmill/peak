function [peak_all, out] = peak_impulse_mimo(A, B, C, order, signed, rank_tol, Tmax)
%PEAK_IMPULSE_MIMO  Find the maximum value of max|C_j x| for each
%input-outpupt impulse response of a MIMO linear system with dynamics 
%       x'=Ax+Bu, y=Cx.
%Input:
%   A,B,C:      Linear System
%   order:      Order of moment relaxation (degree=2*order)
%   signed:     max (Cjx) (signed=1) or |Cjx| (signed=0, default)
%   ranktol:    rank tolerance of moment matrix
%   Tmax:       Maximum time. By default is Infinite (time-independent).
%               Tmax ~= Inf will recover optimal time at optimality, and
%               allow for time-varying safety contours
%Ouputs:
%   peak_val:   Maximum absolute value of impulse response
%   out:        Data structure that contains the peak response information

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

%Parameters of relaxation

if ~iscell(A)
    A = {A};
end

nx = size(A{1}, 1);
ny = size(C, 1);
nu = size(B, 2);
nsys = length(A);
d = order*2;

opt = cell(ny, nu);
peak_val = zeros(ny, nu);

for j = 1:ny
    Cj = C(j, :);        
    for i = 1:nu
        Bi = B(:, i);

        %run peak code on I/O pair (i, j)
        [peak_val_curr, opt_curr] = peak_impulse_siso(A, Bi, Cj, order, signed, rank_tol, Tmax);
        peak_val(j, i) = peak_val_curr;
        opt{j, i} = opt_curr;
    end
end

peak_y = max(peak_val, [], 2);
[peak_all, ind_max] = max(peak_val(:));
[j_max, i_max] = ind2sub([ny, nu], ind_max);
opt_max = opt{j_max, i_max};

%output optimal solution

out = struct;
out.peak_val = peak_val;
out.peak_all = peak_all;
out.peak_y = peak_y;
out.opt = opt;
out.opt_max = opt_max;
out.peak_all = peak_all;
out.ind_max = [j_max, i_max];

end

