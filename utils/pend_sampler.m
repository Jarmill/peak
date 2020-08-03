function [out_sim] = pend_sampler(sampler, Ns, Tmax, nonneg, b, u, odefcn)
%SWITCH_SIM Simulate a system that switches between different dynamics from
%time 0 to Tmax. Switches are modeled by an exponential distribution with
%mean mu. This assumes that all dynamics have the same valid region, will
%need to work on the closed cover further later.
%
% Input:
%   sampler:    A function sampler() that returns a single sampled point on
%               X0 (inital point)
%   Ns:         Number of sampled points to run
%   Tmax:       Maximum time range of simulation
%   nonneg:     Mean time for system switching
%   b:          Friction
%   u:          Control Law
%   odefcn:     Function handle to the ode solver (default ode15s to deal
%               with stiffness)
%
% Output:
%   out:        Data structure with fields
%       t:          time
%       x:          state
%       break_sys:  active system in time break
%       break_time: time breaks for system

if nargin < 3
    Tmax = 10;
end

if nargin < 4
    nonneg = [];
end

if nargin < 5
    b = 0;
end

if nargin < 6
    u = @(x) 0;
end

if nargin < 7
    odefcn = @ode45;
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

    out_sim{i} = pend_sim(x0, Tmax, nonneg, b, u, odefcn);
end

end

