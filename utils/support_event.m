function [event_eval, terminal, direction] = support_event(t, x, supp_eval, Tmin, Tmax)
%SUPPORT_EVENT Event function for switch_sim used when simulating systems
%   return the event function eventFcn. Yes this is tricky


%the event function for support for a system
    %this assumes that t = 0...Tmax-Tmin as with ode45.

    time_supp =  (t + Tmin <= Tmax);

    state_supp = all(supp_eval(reshape(x, [], 1)));

    event_eval = time_supp && state_supp;

    %stop integrating when the system falls outside support
    %
    terminal = 1;
    direction = 1;
end
