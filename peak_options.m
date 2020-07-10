classdef peak_options < handle
    %attributes of peak estimation routine
    %a default set of options
    properties
        
        %% properties of run
        %terminal time   
        T(1,1) double{mustBePositive}  = 1;           
        
        %function dynamics
        f = @(t,x,w) -x; 
        %future feature: have f handle switched systems. f is an array of
        %function handles, and there would be a field X_partition of the
        %switching space        
        
        %objective to minimize
        p function_handle = @(t, x) sum(x.^2);
        
        %% Variables and descriptors
        %variables (array of symbolic variables)
%         t sym = [];     %time
%         x sym = [];     %state
%         w sym = [];     %uncertainty
        var = struct('t', [], 'x', [], 'w', []);
%         P = struct('d', 0, 'e', 1);
        
        %support sets
%         X sym = [];     %all states
%         X0 sym = [];    %initial states
%         XT sym = [];    %final states
%         W sym = [];     %admissable uncertainty
        supp = struct('X', [], 'X0', [], 'XT', [], 'W', []);
        
        %% additional options
        
        %should the measures be time-independent?
        TIME_INDEP(1,1) logical = true; 
        
        %what is the tolerance for identifying a matrix as rank-1?
        rank_tol(1,1) double{mustBePositive} = 1e-3; 
        
        
    end
end