clear
load('sym_attractor_7.mat')
% syms y [2 1];
% vy = out.func.vval(y);
% Lvy = out.func.Lvval(y);
% N = 800;

delta = 0.002;
epsilon = 0.004;

% bnd = @(x, y) 0<x & x<delta & 0<y & y<epsilon;
% vng = @(x) out.func.vval(x) + out.peak_val;
% Lvng = out.func.Lvval;
% 
% initreg = @(x,y) bnd(vng([x; y]), Lvng([x;y]));
% 
% inequalityplot(initreg, [-0.9, 0.9], [-1.5, 1.5])
 
N = 3200;
[XX, YY] = meshgrid(linspace(-0.9, 0.9, N), linspace(-1.5, 1.5, N));  % get 2-D mesh for x and y
z = [XX(:) YY(:)]';
vz = out.func.vval(z);
Lvz = out.func.Lvval(z);

VZ = reshape(vz, N, N);
VZN = VZ + out.peak_val;
LVZ = reshape(Lvz, N, N);

condV = double((0 <= VZN) & (VZN <= delta));
condLV = double((0 <= LVZ) & (LVZ <= epsilon));

cond = condV .* condLV;

[~,] = contourf(XX,YY,cond,[1, 1]); hold on
colormap([ 0.7*[1 1 1] ]);



            scatter(out.x0(1, 1), out.x0(2, 1), 200, 'ob', 'DisplayName', 'Peak Initial', 'LineWidth', 2);        
            scatter(out.xp(1, 1), out.xp(2, 1), 200, '*b', 'DisplayName', 'Peak Achieved', 'LineWidth', 2);        
            scatter(out.x0(1, 2), out.x0(2, 2), 200, 'ob', 'HandleVisibility', 'off', 'LineWidth', 2);        
            scatter(out.xp(1, 2), out.xp(2, 2), 200, '*b', 'HandleVisibility', 'off','LineWidth', 2);        
            legend({'Initial Bounds', 'Peak Initial', 'Peak Achieved'}, 'Location', 'Northwest')