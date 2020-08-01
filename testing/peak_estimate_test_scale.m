%test out arguments of peak_estimate routine
%clear
mset clear
rng(300, 'twister')

%dynamics and support set
mpol('x', 2, 1);

%support
Xsupp = [];

%dynamics
f1 = [x(2); -x(1) + (1/3).* x(1).^3 - x(2)]; %prajna and rantzer
f2 = [-x(1); -x(2)]; %go to origin
f3 = [-0.1*x(1)-x(2); x(1) - 0.1*x(2)]; %decaying spiral
X1 = [];
X2 = [];
X3 = [];
%X3 = [x(1) >= 0];

%f = {f1, f2, f3};
%X = {X1, X2, X3};

f = f1;
X = X1;

if iscell(f)
    nsys = length(f);
else
    nsys = 1;
end


%initial set
C0 = [1.5; 0];
R0 = 0.4;

%C0 = [0.8; 0];
%R0 = 0.1;

X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);

%objective to maximize
objective = -x(2);

p_opt = peak_options;
p_opt.var.x = x;

p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;

%p_opt.Tmax = Tmax_sim;
p_opt.state_init = X0;
p_opt.state_supp = Xsupp;
%p_opt.box = [-6, 6];
p_opt.box = 5;
%p_opt.box = [-3, 4];
%p_opt.box = 6;
p_opt.scale = 0;

p_opt.rank_tol = 4e-3;
p_opt.obj = objective;

order = 4;
out = peak_estimate(p_opt, order);
peak_val = out.peak_val;

p_opt_2 = p_opt;
p_opt_2.scale = 1;
out_scale = peak_estimate(p_opt, order);
peak_val_scale = out_scale.peak_val;