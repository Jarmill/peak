%from the Chesi papers

%switched linear systems

%f1 = @(x) [0 2; -1 -1]*x;
%f2 = @(x) [1 2; -3 -2]*x;
rng(45)
SOLVE = 1;
SAMPLE = 1;
PLOT = 1;

%2d example. Will generalize

A = {[0 2; -1 -1], [1 2; -3 -2]};
nsys = length(A);

%i'm surprised this functional junk works
%A_func is the set of functions x -> Ai x for each Ai in A
A_func = cellfun(@(Ai) (@(t, x) Ai*x), A, 'UniformOutput', false);

B = [1; 1];
%C = [1 3];
% C = [1 0];

nx = size(A{1});
nu = size(B, 2);
ny = size(C, 2);

%will eventually want to compute peak response of this system
%try to evaluate the peak response
if SOLVE
    order = 3;    
    
    [peak_val, opt] = peak_impulse_siso(A, B, C, order);
end

%% Simulations
if SAMPLE
    Tmax = 20;  %Maximum time to plot
    Nsample = 70; %number of trajectories

    %mu = 0.5;   %mean time to switch

    %mu_range = [-2, log10(Tmax)];
    mu_range = [-1.5, log10(Tmax)];
    %some sort of distribution on mu
    mu = 10.^(diff(mu_range)*rand(Nsample, 1) + mu_range(1));
    %mu = 1 * ones(Nsample, 1);



    xtraj = cell(Nsample, 1);
    for i = 1:Nsample    
        xtraj{i} = switch_sim(A_func, B, Tmax, mu(i));
        xcell = num2cell(xtraj{i}.x, 1);
        xtraj{i}.v = opt.pval(xcell{:});
        xtraj{i}.Lv = zeros(length(xtraj{i}.v), nsys);
        for j = 1:nsys
            Lv_curr = opt.Lpval{j}(xcell{:});
            xtraj{i}.Lv(:, j) = Lv_curr;
        end        
        xtraj{i}.cost = (xtraj{i}.x * C').^2;

    end
    
    %evaluate functions along trajectories
end

%% plot the output
if PLOT
    
    FS = 14;
    
    %% States and Contours
    figure(1)
    clf
    hold on
    for i = 1:Nsample
        if i == 1
            plot(xtraj{i}.x(:, 1), xtraj{i}.x(:, 2), 'color', 0.6*[1,1,1], 'DisplayName', 'Switched Trajectory')
        else
            plot(xtraj{i}.x(:, 1), xtraj{i}.x(:, 2), 'color', 0.6*[1,1,1], 'HandleVisibility', 'Off')
        end
    end

    title(['Linear Switched System order = ', num2str(order), ', peak estimate = ', num2str(peak_val, 4)], 'FontSize', FS)
    xlim([-2, 2])
    ylim([-2, 2])
    xlabel('x_1')
    ylabel('x_2')
    scatter(B(1), B(2), 200, 'ok', 'DisplayName', 'Initial Condition', 'LineWidth', 2)
    p_neg = line_range([C, -opt.peak_val], [-2, 2], [-2, 2]);
    p_pos = line_range([C, opt.peak_val],  [-2, 2], [-2, 2]);
    if opt.rankp == 1
        scatter(opt.xp(1), opt.xp(2), 100, '*k', 'DisplayName', 'Peak Value', 'LineWidth', 2)
    end
    if ~isempty(p_neg)
        plot([p_neg{1}(1), p_neg{2}(1)],  [p_neg{1}(2), p_neg{2}(2)], 'r--', 'Linewidth', 3, 'DisplayName', 'Peak Certification')
    end
    if ~isempty(p_pos)
        plot([p_pos{1}(1), p_pos{2}(1)],  [p_pos{1}(2), p_pos{2}(2)], 'r--', 'Linewidth', 3, 'HandleVisibility','off')
    end

    fi2 = fimplicit(opt.p  + opt.obj, [-2, 2, -2, 2], ...
                ':k', 'DisplayName', 'Peak Contour', 'LineWidth', 3, 'MeshDensity', 150);    
    hold off
    axis square
    legend('location', 'southwest')
    
    
    %% Values
    figure(2)
    clf
    tiledlayout(nsys+2, 1);
    
%     subplot(3,1,1)
    nexttile
    hold on       
    for i = 1:Nsample
        if i == 1
            plot(xtraj{i}.t, xtraj{i}.v  + opt.obj, 'c','DisplayName', 'Trajectories')
        else
            plot(xtraj{i}.t, xtraj{i}.v + opt.obj, 'c', 'HandleVisibility','off')
        end  
    end
    
    %plots
    
    plot(xlim, [0, 0], ':k',  'HandleVisibility','off')

    hold off
    title('Safety Function along Trajectories', 'FontSize', FS)
    xlabel('time')
    ylabel('v(x) - \gamma')
    legend('location', 'northwest')
    
%     subplot(3,1,2)
    
    for j = 1:nsys
        nexttile
        hold on
        for i = 1:Nsample
            if i == 1
                plot(xtraj{i}.t, xtraj{i}.Lv(:, j), 'c','DisplayName', 'Trajectories')
            else
                plot(xtraj{i}.t, xtraj{i}.Lv(:, j), 'c', 'HandleVisibility','off')
            end  
        end

        %plots

        plot(xlim, [0, 0], ':k',  'HandleVisibility','off')

        hold off
        title(['Change in Safety Function for System ', num2str(j)], 'FontSize', FS)
        xlabel('time')
        ylabel(['L_{f', num2str(j), '}v(x)'])
        legend('location', 'northwest')
    
    end
%     subplot(3,1,3)  
    nexttile
    hold on
    
    for i = 1:Nsample
        p_comp = -xtraj{i}.cost - xtraj{i}.v;
        
        if i == 1
            plot(xtraj{i}.t, p_comp, 'c','DisplayName', 'Trajectories')
        else
            plot(xtraj{i}.t, p_comp, 'c', 'HandleVisibility','off')
        end  
    end

    %plots
    
    plot(xlim, [0, 0], ':k',  'HandleVisibility','off')
    hold off
    
    
    title('Comparision with Objective Function', 'FontSize', FS)
    xlabel('time')
    ylabel('-cost(x) - v(x,t)')
    legend('location', 'northwest')

    
end