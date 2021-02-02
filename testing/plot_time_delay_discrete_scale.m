options = ddeset('AbsTol', 1e-9, 'RelTol', 1e-7);

T = 1;
tau = 0.1;
k = 10;
sol = ddesd(@(t,y,z) -k*z, @(t,y) t-tau,@(t) -1,[0,T], options);

figure(1)
clf
hold on
plot([-tau sol.x], [-1 sol.y])
plot([-tau, T], [0, 0], ':k')
xlim([-tau, T])
hold off


% function d = curr_delay(t,y)
% 
% d = t - tau;
% end
% 
% function y = history(t)
%     y = -1;
% end
% 
% function dydt = ddex(t,y,Z)
% ylag = Z(:,1);
% % dydt = 0.5*y(1) - ylag(1);
% dydt = -5*ylag(1);
% end