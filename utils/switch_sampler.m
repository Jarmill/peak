function [out_sim] = switch_sampler(dynamics, sampler, Ns, Tmax, mu, Nw, odefcn)
%SWITCH_SIM Simulate a system that switches between different dynamics from
%time 0 to Tmax. Switches are modeled by an exponential distribution with
%mean mu. This assumes that all dynamics have the same valid region, will
%need to work on the closed cover further later.
%
% Turn this into an options datastructure
%
% Input:
%   dynamics:   A struct with fields:
%       f:          System dynamics as a function f: (t, x) -> x' = f(t,x)
%       X:          Support set in space (default everywhere)
%       Tmin:       Minimum time where system is active (default 0)
%       Tmax:       Maximum time where system is active (default Inf)
%   sampler:    A function sampler() that returns a single sampled point on
%               X0 (inital point)
%   Ns:         Number of sampled points to run
%   x0:         Initial condition
%   Nw:         Number of parameters (x0 -> [x0, w0]
%   Tmax:       Maximum time range of simulation
%   mu:         Mean time for system switching
%   odefcn:     Function handle to the ode solver (default ode15s to deal
%               with stiffness)
%
% Output:
%   out:        Data structure with fields
%       t:          time
%       x:          state
%       break_sys:  active system in time break
%       break_time: time breaks for system

if nargin < 4
    Tmax = 10;
end

if nargin < 5
    mu = 1;
end

if nargin < 6
    Nw = 0;
end

if nargin < 7
    odefcn = @ode15s;
end


out_sim = cell(Ns, 1);


for i = 1:Ns
    
    if isa(sampler, 'function_handle')
        x0 = sampler();                
    else
        %numeric
        if size(sampler, 2) == 1
            %single point
            x0 = sampler;
        else
            %array, so Ns = size(sampler, 2)
            x0 = sampler(:, i);
        end
    end

    out_sim{i} = switch_sim(dynamics, x0, Tmax, mu, Nw, odefcn);
end

end

