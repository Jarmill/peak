function [peak_val, xp_val, tp_val, Mp_1] = peak_impulse_siso(A, B, C, order, ranktol, Tmax)
%PEAK_IMPULSE_SISO  Find the maximum value of |Cx| for the impulse response
%of a SISO linear system with dynamics x'=Ax+Bu, y=Cx.
%Input:
%   A,B,C:      Linear System
%   order:      Order of moment relaxation (degree=2*order)
%   ranktol:    rank tolerance of moment matrixof 
%   Tmax:       Maximum time. By default is Infinite (time-independent).
%               Tmax ~= Inf will recover optimal time at optimality, and
%               allow for time-varying safety contours
%Ouputs:
%   peak_val:   Maximum absolute value of impulse response
%   xp_val:     Expectation of optimal point
%   tp_val:     Expectation of optimal time (or Inf if Tmax = Inf)
%   Mp_1:       Top corner of moment matrix

if nargin < 5
    Tmax = Inf;    
end

if nargin < 4
    ranktol = 1e-3;
end

if nargin < 3
    order = 2;
end

d = order*2;
    
outputArg1 = inputArg1;
outputArg2 = inputArg2;
end

