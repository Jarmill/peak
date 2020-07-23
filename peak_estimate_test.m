%test out arguments of peak_estimate routine
clear
%mset clear


%dynamics and support set
%prajna and rantzer flow

%go back to functions
%and/or figure out how to extract monomials and powers from mpol
mpol('x', 2, 1);

%support
Xsupp = [];

%dynamics
f1 = [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];
f2 = [-x(1); -x(2)];
X1 = [];
X2 = [];
f = {f1, f2};
X = {X1, X2};



%dynamics = struct('f', f, 'X', X);
%dynamics = struct('f', {f1, f2}, 'X', {X, X});


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

p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;

%p_opt.Tmax = 10;
p_opt.state_init = X0;
p_opt.state_supp = Xsupp;

p_opt.rank_tol = 4e-3;
p_opt.obj = objective;

order = 3;
out = peak_estimate(p_opt, order);

%now do plots
%out.func has all the ingredients
%will require a for loop

%Things to do:
%write the sampler
%write the plotting code
%further testing

[out_sim] = switch_sim(out.dynamics, C0, 10, 1);

% [ts, ys] = ode45(out.func.fval{1}, [0, 10], C0);
% Nts = length(ts);
% Xs = zeros(Nts, 1);
% for i = 1:Nts
%     Xs(i) = min(out.func.Xval{1}(ys(i, :)'));
% end