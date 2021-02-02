options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-7);

T = 10;
sol = ddesd(@ddex,@curr_delay,@history,[0,T], options);

figure(1)
clf
hold on
plot(sol.x, sol.y)
plot([0, T], [0, 0], ':k')
hold off


function d = curr_delay(t,y)
%DDEX1DELAYS  Delays for using with DDEX1DE.

d = [ 0.2 * t];
% d = t - 1;
end

function y = history(t)
    y = -2;
end

function dydt = ddex(t,y,Z)
%DDEX1DE  Example of delay differential equations for solving with DDE23.
%
%   See also DDE23.

%   Jacek Kierzenka, Lawrence F. Shampine and Skip Thompson
%   Copyright 1984-2014 The MathWorks, Inc.
K = 5;
ylag = Z(:,1);
dydt = K*(-0.5*ylag(1) - y(1));
% dydt = 0.5*y(1) - ylag(1);
% dydt = -ylag(1);
end