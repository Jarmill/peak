function [event_eval, terminal, direction] = support_event(t, x, supp_eval, Tmin, Tmax)
%SUPPORT_EVENT Event function for sampler used when simulating systems
%   return the event function eventFcn. Yes this is tricky


%the event function for support for a system
    %this assumes that t = 0...Tmax-Tmin as with ode45.
if nargin < 5
    time_supp = ones(size(t));
else
    time_supp =  (t + Tmin <= Tmax);
end
    Npt = size(x, 2);
    event_eval = zeros(1, Npt);
    for i = 1:Npt
        xcurr = x(:, i);
        tcurr = t(:, i);
        
        state_supp = all(supp_eval(reshape(xcurr, [], 1)));

        event_eval(i) = time_supp(i) && state_supp;
    end

    %stop integrating when the system falls outside support
    %
    terminal = 1;
    direction = 0;
end
