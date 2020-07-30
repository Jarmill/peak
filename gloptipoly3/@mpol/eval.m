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
nvar = length(var);
%find types of constraints

p_eval= zeros(np, npt);

if isa(pt, 'sym')
    %interface with matlab's symbolic math toolbox
    p_eval = sym(p_eval);
end

%tricky to address relevant variables
varname_cell = cell(nvar, 1);
[varname_cell{:}] = var.var;
%permute input array in terms of gloptipoly's names of variables
[var_name_all, var_perm_all] = sort(cell2mat(varname_cell));

pt2 = pt(var_perm_all, :);
%could probably vectorize this to make it more efficient
for i = 1:np
    curr_pow = p(i).pow;
    curr_coef = p(i).coef;
    curr_var = p(i).var;
    
    %varname_cell = cell(length(curr_var), 1);
    %[varname_cell{:}] = curr_var;
    %permute input array in terms of gloptipoly's names of variables
    %[var_name, var_perm_old] = sort(cell2mat(varname_cell));
    
    if (length(curr_var) == 1) && (curr_var == 0)
        %constant function
        p_eval(i, :) = curr_coef;
    else
        
        %relabel variables
        [curr_var_sort, ~] = sort(curr_var);    
        
        var_perm = find(any(curr_var_sort == var_name_all, 2));
    
        %var_perm = ind_curr_var;
        
        for j = 1:npt
        %p_eval(:, j) = double(subs(p, var, pt(:, j))); 
            pt_trans = reshape(pt2(var_perm, j), 1, length(pt2(var_perm, j)));
            monom = pt_trans.^curr_pow;
            monom_prod = prod(monom, 2);
            p_eval(i, j) = curr_coef' * monom_prod;        
        end
    end
end

