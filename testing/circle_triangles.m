%triagonal approximation to circle

Nth = 15;

th = 2*pi/Nth * (0:(Nth-1));
th_avg = pi/Nth + th;

cth = cos(th);
sth = sin(th);
cth_avg = cos(th_avg);
sth_avg = sin(th_avg);

pts = cell(Nth, 1);

for i = 1:Nth
%for i = 1:1
    pt_curr = zeros(2, 3);
    if i == Nth
        i_next = 1;
    else
        i_next = i + 1;
    end
    
    M = [cth(i), cth(i_next); sth(i) sth(i_next)];
    pt_curr(:, 1:2) = M;
    
%     D = cth(i)*sth(i_next) - cth(i_next)*sth(i);
    
%     pt_intersect = pt_curr(:, 1:2) \ [1; 1];
%     pt_curr(:, 3) = pt_intersect;
    
    
    %pt_curr(:, 3) = [sth(i_next)-sth(i), cth(i_next)-cth(i)]/D;
%     pt_curr(:, 3) = [cth_avg(i), sth_avg(i)];
    
    pt_curr(:, 3) = M' \ [1; 1];

    pts{i} = pt_curr;
end


th_circ = linspace(0, 2*pi, 200);

figure(1)
clf
hold on


for i = 1:Nth
%for i = 1:1
    x_curr = pts{i}(1, :);
    y_curr = pts{i}(2, :);
    
    patch(x_curr, y_curr, 'k')
end    
plot(cos(th_circ), sin(th_circ), 'r', 'linewidth', 3)

xlim([-1.2, 1.2])
ylim([-1.2, 1.2])
axis square
title('Triangular approximation to Circle')