%Example 2.1 from the Fantuzzi/Goluskin paper
%on a circle
mset clear
rng(300, 'twister')

%dynamics and support set
%prajna and rantzer flow

%go back to functions
%and/or figure out how to extract monomials and powers from mpol
mpol('x', 2, 1);

%support
Xsupp = [];

%dynamics
% f = [x(2)*t - 0.1*x(1) - x(1)*x(2);
%     -x(1)*t - x(2) + x(1)^2];

r2 = x'*x;

A = [0.2, 1; 0, -0.4];
J = [0, -1; 1,  0];
f = A*x + J*x*r2;

X = [];

% X = (x(2) >= 0)
%initial set
C0 = [0; 0];
R0 = 0.5;

EQ = 1;
% SYM = 0;

if EQ
    X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 == R0^2);
    sampler = @() sphere_sample(1, 2)'*R0 + C0;
else
    X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);
    sampler = @() ball_sample(1, 2)'*R0 + C0;
end

% if ~SYM
%     sampler0 = sampler;
%     X = [X; V_sym'*x >= 0];
%     sampler = sample_rejector(V_sym, sampler0);
% end


%objective to maximize
% objective = x(1);
%objective = [x(1);  1/sqrt(5) * x(1) + 2/sqrt(5) *x(2)];

objective = r2;

% objective = [x(1); x(2)];

p_opt = peak_options;
p_opt.var.x = x;

p_opt.state_supp = Xsupp;
p_opt.state_init = X0;

p_opt.dynamics = struct;
p_opt.dynamics.f = f;
p_opt.dynamics.X = X;
p_opt.dynamics.discrete = 0;
Tmax_sim = 20;
% p_opt.Tmax = Tmax_sim;

p_opt.box = 2;
p_opt.scale = 0;
%p_opt.R = 6;


p_opt.rank_tol = 3e-4;
p_opt.obj = objective;

order = 7;
out = peak_estimate(p_opt, order);
peak_val = out.peak_val;

%deal with symmetry
M0_mix = out.M0(2:3, 2:3);
Mp_mix = out.Mp(2:3, 2:3);

rank_sym_0 = rank(M0_mix, p_opt.rank_tol);
rank_sym_p = rank(Mp_mix, p_opt.rank_tol);

if rank_sym_0==1 && rank_sym_p==1
    %the solution is optimal with symmetry
    %the symmetry structure of the solution is x <=> -x
    %equivariance: f(-x) = -f(x)
    

    
    %find orbit of best initial point x0
    sym_sign0 = sign(M0_mix(2,1));
    x0_abs = sqrt([M0_mix(1,1); M0_mix(2,2)]);
    x0 = [1; sym_sign0].* x0_abs;
%     x0_neg = -x0_pos;
    
    %find orbit of peak achieved point xp
    sym_signp = sign(Mp_mix(2,1));    
    xp_abs = sqrt([Mp_mix(1,1); Mp_mix(2,2)]);
    xp = [1; sym_signp].* xp_abs;
%     xp_neg = -xp_pos;

    %package into the output
    if ~out.dynamics.time_indep
        out.tp = [1 1] * out.tp;
    end
    out.x0 = [x0 -x0];
    out.xp = [xp -xp];
    out.optimal = 1;
end
    



%% now do plots
%there is something wrong with the definition of v. why?

rng(50, 'twister')
x0 = C0;

mu = 1;

Nsample = 100;

% Nsample = `50;

% Nsample = 20;


out_sim = switch_sampler(out.dynamics, sampler, Nsample, Tmax_sim, mu, 0, @ode45);

if (out.optimal == 1)
    out_sim_peak = switch_sampler(out.dynamics, out.x0, 2, Tmax_sim);
    nplot = nonneg_plot(out, out_sim, out_sim_peak);
    cplot = cost_plot(out, out_sim, out_sim_peak);
    splot = state_plot_2(out, out_sim, out_sim_peak);
else
    nplot = nonneg_plot(out, out_sim);
    cplot = cost_plot(out, out_sim);
    splot = state_plot_2(out, out_sim);
end
