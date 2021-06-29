function [event_eval, terminal, direction] = blank_event(t, x, w)
%BLANK_EVENT an event function that is always true (trivial)

    %stop integrating when the system falls outside support
    %
    event_eval = 1;
    terminal = 1;
    direction = 0;
end