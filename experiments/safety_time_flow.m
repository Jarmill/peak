%safety time by occupation measure
%https://arxiv.org/abs/1903.05311v2

%example for safety margin flow example in fig 6

mset clear;

%% parameters of problem
order = 10;
d = order * 2;

SAFE = 1;

Tmax = 5;
box = 3;

%variables
mpol('x0', 2); mpol('t0');      %initial
mpol('xT', 2); mpol('tT');      %final
mpol('xocc', 2); mpol('tocc');  %occupation
mpol('xs', 2); mpol('ts');      %safe occupation
mpol('xu', 2); mpol('tu');      %unsafe occupation

%measures
m0 = meas([x0; t0]);        mon0 = mmon([x0; t0], d);
mT = meas([xT; tT]);        monT = mmon([xT; tT], d);
mocc = meas([xocc; tocc]);  monocc = mmon([xocc; tocc], d);
ms = meas([xs; ts]);        mons = mmon([xs; ts], d);
mu = meas([xu; tu]);        monu = mmon([xu; tu], d);

%% flow system 
%this flow is time invariant
f = Tmax * [xocc(2); -xocc(1) + (1/3).* xocc(1).^3 - xocc(2)];

%initial set
%C0 = [1.2; 0];
C0 = [1.5; 0];
R0 = 0.4;

%C0 = [0.8; 0];
%R0 = 0.1;

X0 = ((x0(1)-C0(1))^2 + (x0(2)-C0(2))^2 <= R0^2);

%unsafe set

%circle
Cu = [0; -0.5];
%Cu = [2.5; 0];
Ru = 0.5;
c1g = Ru^2 - (xu(1) - Cu(1)).^2 - (xu(2) - Cu(2)).^2;

if SAFE
    theta_c = 5*pi/4;       %p* = -0.1417, beta = [0, 1]
else
    theta_c = 3*pi/4;       %p* = 0.1935, beta = [0.712, 0.287]
end

w_c = [cos(theta_c); sin(theta_c)];
c2g = w_c(1)*(xu(1) - Cu(1)) + w_c(2) * (xu(2) - Cu(2)); 

%unsafe set
g = [c1g >= 0; c2g >= 0];


%% Constraints
supp_con = [X0; t0 == 0;
    xT.^2 <= R^2; tT == 1;
    xocc.^2 <= R^2; tocc*(1 - tocc) >= 0;
    g;  tu*(1 - tu) >= 0;
    xs.^2 <= R^2; ts*(1 - ts) >= 0];

%liouville
Ay = mom(diff(monocc, xocc)*f) + mom(diff(monocc, tocc));
Liou = Ay +  mom(mon0) - mom(monT);

mom_con = [mass(m0) == 1; mom(mons) + mom(monu) == mom(monocc); Liou == 0];

objective = max(mass(mu));


mset('yalmip',true);
mset(sdpsettings('solver', 'mosek'));

P = msdp(objective, ...
    mom_con, supp_con);

[status,obj_rec, m,dual_rec]= msol(P);
time_unsafe = obj_rec * Tmax;