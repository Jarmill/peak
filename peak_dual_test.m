%solve the flow test program in dual form
%hopefully this works


%options and degrees
opts = sdpsettings('solver', 'mosek');
order = 10;
d = 2*order;

%variables in polynomials
x = sdpvar(2, 1);
t = sdpvar(1, 1);

%define polynomial v and parameter gamma
[v, cv] = polynomial([t; x], d);
gamma = sdpvar(1, 1);

%constraint set
R = 5;

%final time
T = 20;

%dynamics
f = T*[x(2); -x(1) + (1/3).* x(1).^3 - x(2)];

Lv = jacobian(v, t) + jacobian(v, x) * f;

% unsafe region level
Ru = 0.4;
Cu = [-1; -1];
p =  Ru^2 - (x(1) - Cu(1))^2 - (x(2)- Cu(2))^2;

%Initial set
R0 = 0.5;
C0 = [1.5; 0];
gX0=  R0^2 - (x(1) - C0(1))^2 + (x(2)- C0(2))^2;
gX = R^2 - x(1)^2 - x(2)^2;
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
%mu0:    v(0, x) >= gamma      x in X0
%mu:     Lv(t,x) >= 0          t in [0,T], x in X
%mup:    p(t,x)  >= v(t,x)     t in [0,T], x in X

%initial conditions
cons = [sos(replace(v, t, 0) - gamma - gX0*sx0), sos(sx0)];
%dynamics
cons = [cons, sos(Lv - gt*stf - gX*sxf), sos(stf), sos(sxf)];
%peak
cons = [cons, sos(p - v - gt*stp - gX*sxp), sos(stp), sos(sxp)];



solvesos(cons, obj, opts, [cv; cx0; ctf; cxf; ctp; cxp; gamma]);

value(obj)