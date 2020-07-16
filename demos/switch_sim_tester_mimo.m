%from the Chesi papers

%switched linear systems

%f1 = @(x) [0 2; -1 -1]*x;
%f2 = @(x) [1 2; -3 -2]*x;
rng(45)
SOLVE = 1;
SAMPLE = 1;
PLOT = 1;

%2d example. Will generalize

%A = {[0 2; -1 -1], [1 2; -3 -2]};
A = {[0 2; -1 -1], [1 2; -3 -2], [-1 0; 1 -0.5]};
nsys = length(A);

%i'm surprised this functional junk works
%A_func is the set of functions x -> Ai x for each Ai in A
A_func = cellfun(@(Ai) (@(t, x) Ai*x), A, 'UniformOutput', false);

%Input
%B = [1; 1];
B = [1 -1; 1 0];
%B = [1 -1 -0.5; 1 0 1];
%B = [1 0.8 1; 1 -0.5 -1];

%Output
%C = [1 3];
C = [1 3; 1 0];
%C = [1 0];
%C = [1 3; -1 0; 0 1; -1 1];
%C = C ./ sqrt(sum(C.^2, 1)); %normalize?

nx = size(A{1}, 1);
nu = size(B, 2);
ny = size(C, 1);

m_plot = 4;

%will eventually want to compute peak response of this system
%try to evaluate the peak response
if SOLVE
    order = 3;
    
    [peak_val, opt] = peak_impulse_mimo(A, B, C, order);

end

%% Simulations
if SAMPLE
    Tmax = 20;  %Maximum time to plot
    
    CORNER_SAMPLE = 20; %number of samples at corners
    Nsammple = nu*CORNER_SAMPLE + 60;
    %Nsample = 70; %number of trajectories
    

    %mu = 0.5;   %mean time to switch

    %mu_range = [-2, log10(Tmax)];
    %mu_range = [-1.5, log10(Tmax)];
    mu_range = [-1.2, log10(Tmax)];
    %mu_range = [-1, log10(Tmax)];
    %some sort of distribution on mu
    mu = 10.^(diff(mu_range)*rand(Nsample, 1) + mu_range(1));
    Bind = randi([1, nu], Nsample, 1);
    
    %mu = 1 * ones(Nsample, 1);



    xtraj = cell(Nsample, 1);
    for i = 1:Nsample       
        if i <= nu * CORNER_SAMPLE
            x0_curr = B(:, mod(i, nu)+1);
        else            
            rand_weights = rand(nu, 1);
            rand_weights = rand_weights/sum(rand_weights);
            x0_curr = B*rand_weights;        
        end
        xtraj{i} = switch_sim(A_func, x0_curr, Tmax, mu(i));
        Nt_curr = length(xtraj{i}.t);
        xcell = num2cell(xtraj{i}.x, 1);
        xtraj{i}.v = zeros(Nt_curr, ny);
        xtraj{i}.Lv = zeros(Nt_curr, ny, nsys);
        xtraj{i}.cost = zeros(Nt_curr, ny);
        for k = 1:ny
            C_curr = C(k, :);
            
            xtraj{i}.v(:, k) = opt.opt{k, Bind(i)}.pval(xcell{:});

            for j = 1:nsys
                
                Lv_curr = opt.opt{k, Bind(i)}.Lpval{j}(xcell{:});
                xtraj{i}.Lv(:, k, j) = Lv_curr;
            end        
            xtraj{i}.cost(:, k) = (xtraj{i}.x * C_curr').^2;
        end
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
    xlim([-m_plot, m_plot])
    ylim([-m_plot, m_plot])
    xlabel('x_1')
    ylabel('x_2')
    %all inputs
    
    
    %initial conditions
    if nu == 2
        plot(B(1, :), B(2, :),'k', 'DisplayName', 'Initial Conditions', 'LineWidth', 2)
    elseif nu > 2
        k = convhull(B');
        plot(B(1, k), B(2, k),'k', 'DisplayName', 'Initial Conditions', 'LineWidth', 2)
    end
    scatter(B(1, :), B(2, :), 200, 'ok', 'DisplayName', 'Pure Input', 'LineWidth', 2)
    
    if opt.opt_max.rankp == 1
        scatter(opt.opt_max.xp(1), opt.opt_max.xp(2), 100, '*k', 'DisplayName', 'Peak Value', 'LineWidth', 2)
    end
    
    
    %output plots
    if ny == 1
        %only one output
        p_neg = line_range([C, -opt.peak_val], [-m_plot, m_plot], [-m_plot, m_plot]);
        p_pos = line_range([C, opt.peak_val],  [-m_plot, m_plot], [-m_plot, m_plot]);
    
        if ~isempty(p_neg)
            plot([p_neg{1}(1), p_neg{2}(1)],  [p_neg{1}(2), p_neg{2}(2)], 'r--', 'Linewidth', 3, 'DisplayName', 'Peak Certification')
        end
        if ~isempty(p_pos)
            plot([p_pos{1}(1), p_pos{2}(1)],  [p_pos{1}(2), p_pos{2}(2)], 'r--', 'Linewidth', 3, 'HandleVisibility','off')
        end
    else
        %more than one output
        %assuming that C has full row rank
        
        %find corners of polytope
        V = qlcon2vert(zeros(nx, 1), [C; -C], [opt.peak_y; opt.peak_y]);
        kV = convhull(V);
        
        Vmax = qlcon2vert(zeros(nx, 1), [C; -C], opt.peak_all * ones(2*ny, 1));
        kVmax = convhull(Vmax);
        
        plot(V(kV, 1), V(kV, 2), 'r:', 'LineWidth', 3, 'DisplayName', 'Polytope Enclosure')
        plot(Vmax(kVmax, 1), Vmax(kVmax, 2), 'r--', 'LineWidth', 3, 'DisplayName', 'Peak Certification')
        
        
    end

    %fi2 = fimplicit(opt.opt_max.p  + opt.opt_max.obj, [-2, 2, -2, 2], ...
    %            ':k', 'DisplayName', 'Peak Contour', 'LineWidth', 3, 'MeshDensity', 150);    
    hold off
    axis square
    legend('location', 'southwest')
    
    
%     %% Values
%     figure(2)
%     clf
%     tiledlayout(nsys+2, 1);
%     
% %     subplot(3,1,1)
%     nexttile
%     hold on       
%     for i = 1:Nsample
%         for u = 1:nu
%             if i == 1
%                 plot(xtraj{i}.t, xtraj{i}.v  + opt.obj, 'c','DisplayName', 'Trajectories')
%             else
%                 plot(xtraj{i}.t, xtraj{i}.v + opt.obj, 'c', 'HandleVisibility','off')
%             end  
%         end
%     end
%     
%     %plots
%     
%     plot(xlim, [0, 0], ':k',  'HandleVisibility','off')
% 
%     hold off
%     title('Safety Function along Trajectories', 'FontSize', FS)
%     xlabel('time')
%     ylabel('v(x) - \gamma')
% %     legend('location', 'northwest')
%     
% %     subplot(3,1,2)
%     
%     for j = 1:nsys
%         nexttile
%         hold on
%         for i = 1:Nsample
%             if i == 1
%                 plot(xtraj{i}.t, xtraj{i}.Lv(:, j), 'c','DisplayName', 'Trajectories')
%             else
%                 plot(xtraj{i}.t, xtraj{i}.Lv(:, j), 'c', 'HandleVisibility','off')
%             end  
%         end
% 
%         %plots
% 
%         plot(xlim, [0, 0], ':k',  'HandleVisibility','off')
% 
%         hold off
%         title(['Change in Safety Function for System ', num2str(j)], 'FontSize', FS)
%         xlabel('time')
%         ylabel(['L_{f', num2str(j), '}v(x)'])
%         legend('location', 'northwest')
%     
%     end
% %     subplot(3,1,3)  
%     nexttile
%     hold on
%     
%     for i = 1:Nsample
%         p_comp = -xtraj{i}.cost - xtraj{i}.v;
%         
%         if i == 1
%             plot(xtraj{i}.t, p_comp, 'c','DisplayName', 'Trajectories')
%         else
%             plot(xtraj{i}.t, p_comp, 'c', 'HandleVisibility','off')
%         end  
%     end
% 
%     %plots
%     
%     plot(xlim, [0, 0], ':k',  'HandleVisibility','off')
%     hold off
%     
%     
%     title('Comparision with Objective Function', 'FontSize', FS)
%     xlabel('time')
%     ylabel('-cost(x) - v(x,t)')
%     legend('location', 'northwest')

    
end