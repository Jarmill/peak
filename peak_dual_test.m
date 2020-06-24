%solve the flow test program in dual form
%hopefully this works


%options and degrees
opts = sdpsettings('solver', 'mosek');
d = 6;

%variables in polynomials
x = sdpvar(2, 1);
t = sdpvar(1, 1);

%define polynomial v and parameter gamma
[v, cv] = polynomial([t; x], d);
gamma = sdpvar(1, 1);

%constraint set
R = 5;

%final time
T = 10;

%dynamics
f = T*[x(2); -x(1) + (1/3).* x(1).^3 - x(2)];

Lv = jacobian(v, t) + jacobian(v, x) * f;

%objective
%maximize, get close to center of unsafe region as possible
px =  0.16 - (x(1)+1)^2 + (x(2)+1)^2;

%Initial set
gX0 = (0.25 - (x(1)-1.5)^2 + x(2)^2);
gX = (R^2 - x'*x);
gt = t*(1-t);


%Putinar Multipliers
[sx0, cx0] = polynomial(x, d-2);
[stf, ctf] = polynomial([t;x], d-2);
[sxf, cxf] = polynomial([t;x], d-2);
[stp, ctp] = polynomial([t;x], d-2);
[sxp, cxp] = polynomial([t;x], d-2);

obj = -gamma;

%dual program:

%cost:   max gamma
%mu0:    gamma - v(0, x) >=0   x in X0
%mu:     Lv(x,t) <= 0          t in [0,T], x in X
%mup:    v(x,t) >= p(x)        t in [0,T], x in X

% %initial conditions
% cons = [sos(gamma - replace(v, t, 0) - gX0*sx0), sos(sx0)];
% %dynamics
% cons = [cons, sos(-Lv - gt*stf - gX*sxf), sos(stf), sos(sxf)];
% %peak
% cons = [cons, sos(v - px - gt*stp - gX*sxf), sos(stp), sos(sxp)];

%initial conditions
cons = [sos(replace(v, t, 0) - gamma - gX0*sx0), sos(sx0)];
%dynamics
cons = [cons, sos(Lv - gt*stf - gX*sxf), sos(stf), sos(sxf)];
%peak
cons = [cons, sos(-px - v - gt*stp - gX*sxf), sos(stp), sos(sxp)];



solvesos(cons, obj, opts, [cv; cx0; ctf; cxf; ctp; cxp], gamma);