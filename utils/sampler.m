function [out_sim] = sampler(dynamics, Np, opts)
%SAMPLER Samples (possibly) multiple trajectories of a system. Optionally
%includes parameter uncertainty.
%
% Turn this into an options datastructure
%
% Input:
%   dynamics:   A struct with fields:
%       f:          System dynamics as a function f: (t, x) -> x' = f(t,x)
%       X:          Support set in space (default everywhere)
%       Tmin:       Minimum time where system is active (default 0)
%       Tmax:       Maximum time where system is active (default Inf)
%   Np:         Number of sampled points to run
%   opts:       A struct of type sampler_options that contains parameters
%               of the sampler
% Output:
%   out:        Data structure with fields
%       t:          time
%       x:          state
%       break_sys:  active system in time break
%       break_time: time breaks for system

out_sim = cell(Np, 1);

if opts.parallel
parfor i = 1:Np
% for i = 1:Np
    if isa(opts.sample.x, 'function_handle')
        x0 = opts.sample.x();      
    else
        %numeric
        %assume that the dimensions of x and w are compatible if arrays are
        %passed in.
        if size(opts.sample.x, 2) == 1
            %single point
            x0 = opts.sample.x;
        else
            %array, so Ns = size(sampler, 2)
            x0 = opts.sample.x(:, i);
        end
    end
    
    if isa(opts.sample.w, 'function_handle')
        w0 = opts.sample.w();
    else        
        if size(opts.sample.w, 2) == 1
            %single point            
            w0 = opts.sample.w;
        else
            %array, so Ns = size(sampler, 2)            
            w0 = opts.sample.w(:, i);
        end
    end

    %do discrete later
    
    
    %old code
    if isfield(dynamics, 'discrete') && dynamics.discrete
        out_sim{i} = sampler_discrete(dynamics, x0, w0, opts);
    else
        out_sim{i} = sampler_cont(dynamics, x0, w0, opts);
    end
end
else
    for i = 1:Np
% for i = 1:Np
    if isa(opts.sample.x, 'function_handle')
        x0 = opts.sample.x();      
    else
        %numeric
        %assume that the dimensions of x and w are compatible if arrays are
        %passed in.
        if size(opts.sample.x, 2) == 1
            %single point
            x0 = opts.sample.x;
        else
            %array, so Ns = size(sampler, 2)
            x0 = opts.sample.x(:, i);
        end
    end
    
    if isa(opts.sample.w, 'function_handle')
        w0 = opts.sample.w();
    else        
        if size(opts.sample.w, 2) == 1
            %single point            
            w0 = opts.sample.w;
        else
            %array, so Ns = size(sampler, 2)            
            w0 = opts.sample.w(:, i);
        end
    end

    %do discrete later
    
    
    %old code
    if dynamics.discrete
        out_sim{i} = sampler_discrete(dynamics, x0, w0, opts);
    else
        out_sim{i} = sampler_cont(dynamics, x0, w0, opts);
    end
end
end
end

