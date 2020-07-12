%sys = ss([0 1; -0.5 -1], [0 1; 1 1], [1 0;  0 1], 0);
rng(42);

%522, 643 done

% nx = 6;
% ny = 4;
% nu = 3;
% sys = rss(nx, ny, nu);


%Jet Aircraft from MIMO Matlab tutorial
A = [-0.0558   -0.9968    0.0802    0.0415
      0.5980   -0.1150   -0.0318         0
     -3.0500    0.3880   -0.4650         0
           0    0.0805    1.0000         0];

B = [ 0.0073         0
     -0.4750    0.0077
      0.1530    0.1430
           0         0];

C = [0     1     0     0
     0     0     0     1];

D = [0     0
     0     0];

states = {'beta' 'yaw' 'roll' 'phi'};
inputs = {'rudder' 'aileron'};
outputs = {'yaw rate' 'bank angle'};

sys = ss(A,B,C,D,'statename',states,...
'inputname',inputs,...
'outputname',outputs);

nx = size(sys.A, 1);
ny = size(sys.C, 1);
nu = size(sys.B, 2);


SOLVE = 0;
PLOT_CONTOUR = 0;
PLOT_VALUES  = 1;
PLOT = PLOT_CONTOUR || PLOT_VALUES;

Tmax_plot = 15;
trange = [0, Tmax_plot];

%assume all entries of B, C are nontrivial
if SOLVE
    order = 2;
    
    [peak_val, out] = peak_impulse_mimo(A, B, C, order);        
end

if PLOT
    [xtraj.y, xtraj.t, xtraj.x] = impulse(sys, Tmax_plot);
    
    
    %Impulse response plot
    if PLOT_CONTOUR
    figure(1)
    clf
    tiledlayout(nu, ny);
    
    
    for i = 1:nu
        for j = 1:ny
            nexttile
            hold on
            plot(xtraj.t, xtraj.y(:, j, i));
            plot(trange, [0, 0], ':k', 'HandleVisibility', 'Off')
            plot(trange, out.peak_all*[1,1], 'r--', 'Linewidth', 3, 'DisplayName', 'Peak Certification')
            plot(trange, -out.peak_all*[1,1], 'r--', 'Linewidth', 3, 'HandleVisibility', 'off')
            plot(trange, out.peak_val(j, i)*[1,1], 'r:', 'Linewidth', 2, 'DisplayName', 'Current Peak')
            plot(trange, -out.peak_val(j, i)*[1,1], 'r:', 'Linewidth', 2, 'HandleVisibility', 'off')
            hold off
            
            
            if isempty(sys.OutputName{1,1}) || isempty(sys.InputName{1,1})
                title_str = ['Input ', num2str(i), ' to Output ', num2str(j), ', peak = ', num2str(out.peak_val(j, i), 4)];
            else
                title_str = ['Input ', sys.InputName{i}, ' to Output ', sys.OutputName{j}, ', peak = ', num2str(out.peak_val(j, i), 4)];
            end
            
            if [j, i] == out.ind_max
                title_str = [title_str, ' (optimal)'];
            end
            
            title(title_str)
            xlabel('time')
            ylabel('output')
        end
    end
    end
    %nonnegative values
            title(title_str)
    
    xtraj.v = zeros(length(xtraj.t), ny, nu);
    xtraj.Lv = zeros(size(xtraj.v));
    xtraj.cost = zeros(size(xtraj.v));
    for i = 1:nu
        
        xcurr = xtraj.x(:, :, i);       
        xcell = num2cell(xcurr, 1);
        for j = 1:ny
            xtraj.v(:, j, i) = out.opt{j, i}.pval(xcell{:}) + out.opt{j, i}.obj;
            xtraj.Lv(:, j, i) = out.opt{j, i}.Lpval(xcell{:});
            xtraj.cost(:, j, i) = out.opt{j, i}.cost(xcurr);
        end       
    end
    
    if PLOT_VALUES
        figure(2)
        clf
        subplot(3, 1, 1)
        hold on
        for i = 1:nu
            for j = 1:ny
                if [j, i] == out.ind_max
                    plot(xtraj.t, xtraj.v(:, j, i) + out.opt{j, i}.obj, 'b', 'LineWidth', 2, 'DisplayName', 'Peak Traj.')           
                else
                    plot(xtraj.t, xtraj.v(:, j, i) + out.opt{j, i}.obj, 'c', 'HandleVisibility', 'off')
                end                
            end
        end
        hold off
        
        title('Safety Function along Trajectories')
        xlabel('time')
        ylabel('v(x) - \gamma')
        
        subplot(3, 1, 2)
        hold on
        for i = 1:nu
            for j = 1:ny
                if [j, i] == out.ind_max
                    plot(xtraj.t, xtraj.Lv(:,j ,i), 'b', 'LineWidth', 2, 'DisplayName', 'Peak Traj.')           
                else
                    plot(xtraj.t, xtraj.Lv(:,j ,i), 'c', 'HandleVisibility', 'off')
                end                
            end
        end
        hold off
        title('Change in Safety Function along Trajectories')
        xlabel('time')
        ylabel('L_f v(x)')
        
        
        subplot(3, 1, 3)
        hold on
        for i = 1:nu
            for j = 1:ny
                if [j, i] == out.ind_max
                    plot(xtraj.t, xtraj.cost(:,j ,i), 'b', 'LineWidth', 2, 'DisplayName', 'Peak Traj.')           
                else
                    plot(xtraj.t, xtraj.cost(:,j ,i), 'c', 'HandleVisibility', 'off')
                end                
            end
        end
        hold off
        title('Comparision with Objective Function')
        xlabel('time')
        ylabel('-cost(x) - v(x,t)')
    
        
    end
        
end