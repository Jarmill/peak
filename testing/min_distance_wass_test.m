%use a moment relaxation to find the minimum distance between two sets

mset clear
mset('yalmip',true);
mset(sdpsettings('solver', 'mosek'));

order = 2;
d = 2*order;

box_size = 3;


% mpol('x', 2)    %x base
mpol('xw', 2)   %x wasserstein
mpol('xu', 2)   %x unsafe

eta = meas([xw; xu]);

% v = mmon([xw; xu] ,d);

c = mom(sum((xw-xu).^2));

% u_loc = [-1; -1];
% x_loc = [2; -1];
x_loc = [0.5; 2];
% supp_con =  [xw==x_loc; xu==u_loc];
% supp_con = [xw==x_loc; sum(xu.^2)<=1];
supp_con = [xw==x_loc; sum(xu.^2)<=1; xu(2) <= 0];
mom_con = (mass(eta)==1);
objective = min(c);

%pose LMI
P = msdp(objective, ...
    mom_con, supp_con);
    
%solve LMI moment problem
[status, obj] = msol(P);    

Meta = double(mmat(eta));
disp(double(mom(xu)))