%safety margins vs. distance of half-circle set

R = 1;
% x_in = [-0.5; -0.5];
% dist_out = half_circ_dist(x_in, R)

%% plot contours
Ntheta =  75;
clist = 0.1:0.1:2;
% x = dist_contour(Ntheta, R, 0.5);
figure(1)

w = [1, 1];
x_safe = safety_contour(2*Ntheta, R, -0.5, w);
x0 = dist_contour(Ntheta, R, 0);

clf
% plot(x_safe(1, :), x_safe(2, :))
ax1 = subplot(1, 2, 1);
hold on
patch(x0(1, :), x0(2, :), 'r', 'Linewidth', 3, 'EdgeColor', 'none')
title('Contours of L2 Distance')
xlabel('x')
ylabel('y')

ax2  = subplot(1, 2, 2);
hold on
patch(x0(1, :), x0(2, :), 'r', 'Linewidth', 3, 'EdgeColor', 'none')
title('Contours of Safety Margin')
xlabel('x')
ylabel('y')

for i = 1:length(clist)
    subplot(1 ,2, 1);
    x = dist_contour(Ntheta, R, clist(i));
    plot(x(1, :), x(2, :), 'r')
    
    subplot(1 ,2, 2);
    x = safety_contour(Ntheta, R, -clist(i));
    plot(x(1, :), x(2, :), 'r')
end
% title('Contours of L2 Distance')

linkaxes([ax1; ax2])




function dist_out = half_circ_dist(x_in, R)
    %return the L2 distance  between the point x_in and the half circle
    %||x_in||^2 <= R^2 intersect x_in(2) <= 0.
    
    if x_in(2) >= 0
        %flat region
        if x_in(1) < -R
            dist_out = hypot(x_in(1)+R, x_in(2))^2;
        elseif x_in(1) > R
            dist_out = hypot(x_in(1)-R, x_in(2))^2;
        else
            dist_out = x_in(2)^2;
        end
    else
        %circle region
        dist_out = max(norm(x_in, 2)-R, 0);
    end

end

function x_dist = dist_contour(Ntheta, R, c)
    %compute a contour at distance c away from the half-circle with N_theta
    %sample points


    theta_q1 = linspace(0, pi/2, Ntheta);
    theta_q2 = linspace(pi/2, pi, Ntheta);
    theta_q34 = linspace(pi, 2*pi, 2*Ntheta);

    %contour level
    

    x_right = [c*cos(theta_q1)+R; c*sin(theta_q1)];
    x_left = [c*cos(theta_q2)-R; c*sin(theta_q2)];
    x_bottom = [(c+R)*cos(theta_q34); (c+R)*sin(theta_q34)];
    x_dist = [x_right, x_left, x_bottom];
end

function x_safe = safety_contour(Ntheta, R, c, w)
%add constraint weighting later
if nargin < 4
    w = [1; 1];
end

theta_c = 3*pi/2;

radius_out = sqrt(R^2 - c);
if (-c) >= radius_out
    %line doesn't matter anymore
    theta = linspace(0, 2*pi, 3*Ntheta);
    Xmargin = (radius_out)*[cos(theta); sin(theta)];
    x_safe = [Xmargin];

else
phi = asin(-c/radius_out);

th_new_top = phi + theta_c + pi/2; 
th_new_bot = pi - phi + theta_c + pi/2;
 
if th_new_top > pi
    th_new_top = th_new_top - 2*pi;
    th_new_bot = th_new_bot - 2*pi;
end

theta_new_range = linspace(th_new_bot, th_new_top + 2*pi, Ntheta);
circ_margin = [cos(theta_new_range); sin(theta_new_range)];
Xmargin = circ_margin*radius_out;
Xmargin_corner = circ_margin([1,2],[1,Ntheta])*radius_out;

x_safe = [Xmargin, Xmargin_corner];
end

end