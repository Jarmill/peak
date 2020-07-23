function [p_eval] = eval(p, var, pt)
% @MPOL/EVAL Evaluate the polynomial p at point pt
% 
% EVAL(P, VAR, PTS) given the polynomial P of type mpol,
%   evaluate P in variables VAR at point PT of type double
%
%
%
% Does not require @mpol/assign to fix point measures
% useful to keep the variables free
%
% Assume P is a scalar or column vector (edit later)
% Revise this to handle matrix valued p in future
%
% J. Miller, 22 July 2020

npt = size(pt, 2);
np = length(p);
%find types of constraints

p_eval= zeros(np, npt);

for j = 1:npt
    p_eval(:, j) = double(subs(p, var, pt(:, j)));    
end

