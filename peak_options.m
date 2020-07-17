classdef peak_options < handle
    %attributes of peak estimation routine
    %a default set of options
    properties
        
        %% properties of run
        %terminal time   
        Tmax(1,1) double{mustBePositive}  = Inf;           
        
        %function dynamics
        %f: dynamics
        %X: space on which dynamics are valid (arbitrary switching)
        dynamics = struct('f', [], 'X', {})
        
        
        %objective to minimize
        %could be a function (single objective)
        %or a cell of functions (minimum of entries)        
        obj = [];
        
        %% Variables and descriptors
        %variables (array of symbolic variables)
%         t sym = [];     %time
%         x sym = [];     %state
%         w sym = [];     %uncertainty
        var = struct('t', [], 'x', [], 'w', []);
        
        %support sets
%         X sym = [];     %all states
%         X0 sym = [];    %initial states
%         XT sym = [];    %final states (implement time breaks later)
%         W sym = [];     %admissable uncertainty
        %supp = struct('X', [], 'X0', [], 'XT', [], 'W', []);        
        
        
        %(X, T): Xt is the support of trajectories at time Tt
        %W:      Plausible Uncertainties (default to 0)
        state_fix = struct('X', [], 'T', 0)
        
        param = [];        
        
        %% additional options
        
        %should the measures be time-independent?
        %TIME_INDEP(1,1) logical = true; 
        
        %what is the tolerance for identifying a matrix as rank-1?
        rank_tol(1,1) double{mustBePositive} = 1e-3; 
        
        
    end
end