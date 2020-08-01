function [Xfunc] = constraint_func(X,vars, eq_tol)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if nargin < 3
    eq_tol = 1e-8;
end

X_ineq = polyval_func(X.ineq, vars);
X_eq = polyval_func(X.eq, vars);

eq_tol = 1e-8;
Xfunc= @(vart) all([X_ineq(vart) >= 0; abs(X_eq(vart)) <= eq_tol]);
    
end

