%find the minimum distance between trajectories of the 'flow' system and a
%half-circle unsafe set

mset clear
mset('yalmip',true);
mset(sdpsettings('solver', 'mosek'));

order = 2;
d = 2*order;

%% problem parameters

Tmax = 5;

mpol('x', 2, 1);
f = Tmax * [x(2); -x(1) + (1/3).* x(1).^3 - x(2) ];

%initial set
%C0 = [1.2; 0];
C0 = [1.5; 0];
R0 = 0.4;
X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);

%half-circle set
theta_c = 5*pi/4;       %p* = -0.1417, beta = [0, 1]
w_c = [cos(theta_c); sin(theta_c)];
c2f = w_c(1)*(x(1) - Cu(1)) + w_c(2) * (x(2) - Cu(2)); 

BOX = 3;

%% set up measures
mpol('t0', 1, 1);
mpol('x0', 2, 1);
mu0 = meas([t0; x0]);

mpol('t_occ', 1, 1);
mpol('x_occ', 2, 1);
mu_occ = meas([t_occ; x_occ]);

mpol('tp', 1, 1);
mpol('xp', 2, 1);
mup = meas([tp; xp]);
%wasserstein
mpol('xw', 2, 1)
mpol('xu', 2, 1)
eta = meas([xw; xu]);


%% support constraints
supp_con = [t0 == 0; t_occ*(1-t_occ)>=0; tp*(1-tp) >=0;
    