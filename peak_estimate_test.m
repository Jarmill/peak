%test out arguments of peak_estimate routine




clear
mset clear



%define variables
%syms t
%syms('x', [2, 1]);

%handle symbolic toolbox later
%sym2mpol
mpol('x', 2, 1);

%dynamics and support set
%prajna and rantzer flow

%go back to functions
%and/or figure out how to extract monomials and powers from mpol
f = [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];
%f = @(x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];

R = 5;
X = (sum(x.^2) <= R^2);

dynamics.f = f;
dynamics.X = X;

%initial set
C0 = [1.5; 0];
R0 = 0.4;

X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);

%objective to maximize
objective = -x(2);
%objective = -x(2) - x(1);
%
p_opt = peak_options;
p_opt.var.x = x;

p_opt.dynamics = dynamics;
%p_opt.dynamics.f = f;
%p_opt.dynamics.X = X;


p_opt.state_init = X0;
p_opt.state_supp = X;
%p_opt.state_fix.X = X0;
%p_opt.state_fix.X = {X0, X0};
%p_opt.state_fix.T = [0, 1];
p_opt.rank_tol = 4e-3;
p_opt.obj = objective;

order = 4;
out = peak_estimate(p_opt, order);

%now do plots