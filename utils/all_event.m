function [event_eval, terminal, direction] = all_event(t, x)
%SUPPORT_ALL Event function for switch_sim used when simulating systems
%   return the event function eventFcn. Yes this is tricky
%is always true.

    %stop integrating when the system falls outside support
    %
    event_eval = 1;
    terminal = 1;
    direction = 0;
end
