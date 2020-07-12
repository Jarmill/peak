%from the Chesi papers

%switched linear systems

%f1 = @(x) [0 2; -1 -1]*x;
%f2 = @(x) [1 2; -3 -2]*x;
rng(45)

A = {[0 2; -1 -1], [1 2; -3 -2]};

%i'm surprised this functional junk works
%A_func is the set of functions x -> Ai x for each Ai in A
A_func = cellfun(@(Ai) (@(t, x) Ai*x), A, 'UniformOutput', false);

B = [1; 1];
C = [1 3];

%will eventually want to compute peak response of this system

Tmax = 20;  %Maximum time to plot
Nsample = 50; %number of trajectories

%mu = 0.5;   %mean time to switch

mu_range = [-2, log10(Tmax)];
%some sort of distribution on mu
mu = 10.^(diff(mu_range)*rand(Nsample, 1) + mu_range(1));
%mu = 1 * ones(Nsample, 1);



out = cell(Nsample, 1);
for i = 1:Nsample    
    out{i} = switch_sim(A_func, B, Tmax, mu(i));
end

%plot the output

figure(1)
clf
hold on
for i = 1:Nsample
    if i == 1
        plot(out{i}.x(:, 1), out{i}.x(:, 2), 'k', 'DisplayName', 'Switched Trajectory')
    else
        plot(out{i}.x(:, 1), out{i}.x(:, 2), 'k', 'HandleVisibility', 'Off')
    end
end

title('Switched System Plots')
xlim([-2, 2])
ylim([-2, 2])
xlabel('x_1')
ylabel('x_2')
scatter(B(1), B(2), 200, 'ok', 'DisplayName', 'Initial Condition', 'LineWidth', 2)
hold off
legend('location', 'northwest')