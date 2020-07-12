function [peak_all, out] = peak_impulse_mimo(A, B, C, order, rank_tol, Tmax)
%PEAK_IMPULSE_MIMO  Find the maximum value of max|C_j x| for each
%input-outpupt impulse response of a MIMO linear system with dynamics 
%       x'=Ax+Bu, y=Cx.
%Input:
%   A,B,C:      Linear System
%   order:      Order of moment relaxation (degree=2*order)
%   ranktol:    rank tolerance of moment matrixof 
%   Tmax:       Maximum time. By default is Infinite (time-independent).
%               Tmax ~= Inf will recover optimal time at optimality, and
%               allow for time-varying safety contours
%Ouputs:
%   peak_val:   Maximum absolute value of impulse response
%   out:        Data structure that contains the peak response information

if nargin < 6
    Tmax = Inf;    
end

if nargin < 5
    rank_tol = 1e-3;
end

if nargin < 4
    order = 2;
end

%Parameters of relaxation

nx = size(A, 1);
ny = size(C, 1);
nu = size(B, 2);
d = order*2;

opt = cell(ny, nu);
peak_val = zeros(ny, nu);

for j = 1:ny
    Cj = C(j, :);        
    for i = 1:nu
        Bi = B(:, i);

        %run peak code on I/O pair (i, j)
        [peak_val_curr, opt_curr] = peak_impulse_siso(A, Bi, Cj, order);
        peak_val(j, i) = peak_val_curr;
        opt{j, i} = opt_curr;
    end
end

[peak_all, ind_max] = max(peak_val(:));
[j_max, i_max] = ind2sub([ny, nu], ind_max);
opt_max = opt{j_max, i_max};

%output optimal solution

out = struct;
out.peak_val = peak_val;
out.opt = opt;
out.opt_max = opt_max;
out.peak_all = peak_all;
out.ind_max = [j_max, i_max];

end

