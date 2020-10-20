classdef sampler_options < handle
    %attributes of sampling code
    %sample trajectories from (possibly uncertain) systems
    %a default set of options
    properties
        
        %% properties of run
        %terminal time of sampler 
        Tmax(1,1) double{mustBePositive}  = 10;           
        
        %expected time to switch systems (continuous time)
        %exponential distribution with parameter mu
        mu = 1;
        
        %number of variables
        Nw = 0; %time-independent parameters
        Nd = 0; %time-dependent general uncertainty
        Nb = 0; %time-dependent box uncertainty
        
        %function handle
        odefcn = @ode15s;
        
        
        %samplers (default, return empty array)
        sample = struct('x', @() [], ...
                        'w', @() [], ...
                        'd', @() []);
        
    end
end