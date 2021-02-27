function newp = subs_vars(p, old, new)
% @MPOL/SUBS_VAR substitutes variables in a given polynomial to another set 
% of variables. Based on @MPOL/SUBS
%
% SUBS(P, OLD, NEW) given the polynomial P of type mpol, substitutes
%   the variable (or vector of variables) OLD with the polynomial (or
%   vector of polynomials) NEW

% J. Miller, 27 May 2021

%% Input Handling

if length(old)~=length(new)
    error('the second and third argument must have equal size');
end

%% Substitution



[nrows, ncols] = size(p);
newp = p;


if ~isnumeric(p)
    v_old = indvar(old);
    v_new = indvar(new);

    for rind = 1:nrows
        for cind = 1:ncols
    %     p.var
            v_p   = p(rind, cind).var;
            %credit to Guillaume https://www.mathworks.com/matlabcentral/answers/401395-function-changem-or-substitute-values-of-a-matrix#answer_320811
            %for the 'changem' mockup. Ismember may be slow, will need to
            %see.
            [toreplace, bywhat] = ismember(v_p, v_old);
            v_p(toreplace) = v_new(bywhat(toreplace));
            newp(rind, cind).var = v_p;
        end
    end
end