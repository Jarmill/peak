function newsc = subs_vars(sc, old, new)
% @SUPCON/SUBS_VARS substitutes variables in a given support constraint
% Calls @MPOL/SUBS_VARS, faster and with less flexibility/protection than
% @SUPCON/SUBS for the case of a change in basis variables. 
%
% Example: x^2 <= 1 --> y^2 <= 1
%
% SUBS_VARS(SC, OLD, NEW) given the support constraint SC of type supcon, substitutes
%    the vector of variables OLD with the vector of polynomials NEW
%
% J. Miller, 27, Feb 2021
%
% Based on @SUPCON/SUBS by 
% C. Savorgnan, 14 May 2008


newsc=[];

for index=1:length(sc)
    left = subs_vars(sc(index).left, old, new);
    right = subs_vars(sc(index).right, old, new);
    type = sc(index).type;
    newsc = [newsc; supcon(left, right, type)];
end
