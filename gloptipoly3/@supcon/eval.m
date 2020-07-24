function [sc_eval, ineq, eq] = eval(sc, var, pt, tol)
% @SUPCON/EVAL Evaluate the support constraint at point pt
% 
% EVAL(SC, VAR, PTS) given the support constraint of SC of type supcon,
%   evaluate SC in variables VAR at point PT of type double
%
% INEQ: values of inequality constraints at PT (assumed >= 0)
% EQ:   values of equality constraints at PT
% tol: tolerance
%
% J. Miller, 22 July 2020

if nargin < 4
    tol = 1e-10;
end

npt = size(pt, 2);

%find types of constraints
nineq = 0;
neq = 0;
for i = 1:length(sc)
    if strcmp(sc(i).type, 'eq')
        neq = neq + 1;
    else
        nineq = nineq + 1;
    end
end

ineq = zeros(nineq, npt);
eq = zeros(neq, npt);

sc_eval = zeros(length(sc), npt);


for j = 1:npt
    %sc_subs = subs(sc, var, pt(:, j));
    dl = eval(sc.left, var, pt(:, j));
    dr = eval(sc.right, var, pt(:, j));
    for i = 1:length(sc)    
        type = sc(i).type;

        %dl = double(sc_subs(i).left);
        %dr = double(sc_subs(i).right);
        if strcmp(type, 'le')
            sc_eval(i, j) = logical(dl(i) <= dr(i));
            ineq(i, j) = dr(i) - dl(i);
        elseif strcmp(type, 'ge')
            sc_eval(i, j) = logical(dl(i) >= dr(i));
            ineq(i, j) = dl(i) - dr(i);
        else %equal
            sc_eval(i, j) = logical(abs(dl(i)-dr(i)) <= tol);
            eq(i, j) = dl(i) - dr(i);
        end
    end
end
end

