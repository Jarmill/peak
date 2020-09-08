function [out_sim] = switch_sim_discrete(dynamics, x0, Tmax, Nw)
%SWITCH_SIM Simulate a discrete system that switches between different 
%dynamics from integer time 0 to Tmax. Valid regions may be different 
%between subsystems. A random subsystem is chosen for each time index.

if nargin < 4
    Nw = 0;
end

%gather information about the system

%time range of states
if iscell(dynamics.f)
    nsys = length(dynamics.f);    
else
    nsys = 1;
    dynamics.f = {dynamics.f};    
    dynamics.event = @(t,x) deal(1, 1, 0);
end

%support in space and time
nx = length(x0) - Nw;


T = (0:Tmax);
X = zeros(Tmax + 1, nx);
X(1, :) = x0;
x_curr = x0;
system_choice = [];

for k = 2:(Tmax+1)
    
    %find which systems are valid
    possible_sys = [];
    for i = 1:nsys        
        if Nw > 0
            event_value = dynamics.event{i}(k, x_curr, w0);
        else
            event_value = dynamics.event{i}(k, x_curr);
        end
        if event_value == 1
            possible_sys = [possible_sys; i];
        end
    end
    N_possible = length(possible_sys);
    
    if N_possible == 0
        break
    end
    
    curr_sys_ind = randi([1, N_possible], 1, 1);
    
    %choose which system to use. curr_event not needed, validity evaluated
    %on next state iteration
    curr_sys = possible_sys(curr_sys_ind);
    system_choice(k-1) = curr_sys;
    if Nw > 0
        %includes uncertain fixed parameters
        curr_f = @(t, x) dynamics.f{curr_sys}(t, x, w0);        
    else
        %no uncertain fixed parameters
        curr_f = dynamics.f{curr_sys};        
    end
    
    %take the step and store results
    x_next = curr_f(k, x_curr);
    X(k, :) = x_next;
    
    x_curr = x_next;    
end

%output data structure
out_sim = struct;
%system trajectories
out_sim.t = T;
out_sim.x = X;
if Nw > 0
    out_sim.w = w0;
else
    out_sim.w = [];
end

%switching and statistics

out_sim.system = system_choice;
out_sim.Tmax = Tmax;
out_sim.cost = dynamics.cost(X');
%out_sim.nonneg = dynamics.nonneg(X');

if isfield(dynamics, 'nonneg')    
    if Nw > 0        
        out_sim.nonneg = dynamics.nonneg(X', w0);        
    else
        out_sim.nonneg = dynamics.nonneg(X');        
    end
end

end

