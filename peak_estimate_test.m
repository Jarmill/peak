%test out arguments of peak_estimate routine
%clear
%mset clear
rng(300, 'twister')

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
f3 = [-0.1*x(1)-x(2); x(1) - 0.1*x(2)];
%X1 = [x(2) <= 2];
X1 = [];
X2 = [];
%X3 = [];
X3 = [x(1) >= 0];
f = {f1, f2, f3};
X = {X1, X2, X3};

if iscell(f)
    nsys = length(f);
else
    nsys = 1;
end

% f = {f1, f3};
% X = {X1, X3};
%

%f = {f3, f2};
%X = {X3, X2};

% f = f1;
% X = X1;

%f = f3;
%X = X3;


%dynamics = struct('f', f, 'X', X);
%dynamics = struct('f', {f1, f2}, 'X', {X, X});


%initial set
%C0 = [1.2; 0];
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
p_opt.R = 6;


p_opt.rank_tol = 4e-3;
p_opt.obj = objective;

order = 4;
out = peak_estimate(p_opt, order);
peak_val = out.peak_val;


%% now do plots
%there is something wrong with the definition of v. why?

rng(50, 'twister')
x0 = C0;
Tmax_sim = 15;
mu = 1;

Nsample = 100;
sampler = @() circle_sample(1)'*R0 + C0;



out_sim = switch_sampler(out.dynamics, sampler, Nsample, Tmax_sim);

if (out.optimal == 1) && (nsys == 1)
    out_sim_peak = switch_sampler(out.dynamics, out.x0, 1, Tmax_sim);
    nplot = nonneg_plot(out_sim, out_sim_peak);
else
    nplot = nonneg_plot(out_sim);
end


% if iscell(f)
%     nsys = length(f);
% else
%    nsys = 1;
%end


% figure
% clf
% for i = 1:Nsample
%     t_curr = out_sim{i}.t;
%     nonneg_curr = out_sim{i}.nonneg;
%     for k = 1:(nsys+2)
%         subplot(nsys+2, 1, k)
%         hold on
%         plot(t_curr, nonneg_curr(k, :), 'c')
%         
%         if i == 1
%             if k == 1
%                 title('Value Function')
%             elseif k == nsys+2
%                 title('Cost Comparision')
%             else
%                 title(['Lie Derivative across System ', num2str(k-1)])
%             end
%         elseif i == Nsample
%             plot([0, Tmax_sim], [0, 0], ':k')
%         end
%     end       
% end
% 
% if out.optimal    
%     for k = 1:(nsys+2)
%         subplot(nsys+2, 1, k)
%         hold on
%         plot(out_sim_peak{1}.t, out_sim_peak{1}.nonneg(k, :), 'b', 'LineWidth', 3)       
%     end       
% end



% [out_sim] = switch_sim(out.dynamics, x0, Tmax_sim, mu, @ode45);
% 
% figure(1)
% clf
% hold on
% plot(out_sim.x(:, 1), out_sim.x(:, 2))
% scatter(x0(1), x0(2), [], 200, 'ko')
% 
% %implicit plot
% syms tc xc yc;
% if out.dynamics.time_indep
%     vv = conj(monolistYalToGlop([xc; yc], 2*order));
% end
% vsym = out.func.dual_rec'*vv;
% 
% fi2 = fimplicit(vsym  - out.peak_val, [-0.5, 1.5, -2, 2], ...
%             ':k', 'DisplayName', 'Peak Contour', 'LineWidth', 3, 'MeshDensity', 150);    
% 
% 
% figure(2)
% clf
% 
% nn = out_sim.nonneg;
% Nnn = size(out_sim.nonneg, 1);
% %tiledlayout(Nnn + 2, 1);
% for i = 1:Nnn
%     subplot(Nnn, 1, i)
%     hold on
%     plot(out_sim.t, out_sim.nonneg(i, :));
%     plot(out_sim.t([1, end]), [0, 0], ':k')
%     hold off
%     %nexttile;
% end