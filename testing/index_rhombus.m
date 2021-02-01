r = 6;

% i = 0:r;
% j = -r:r;

xb = []; %boundary
xi = []; %interior

for i = 0:r
    for j = -r:r
        if (1 <= (i+j)) && ((i+j) <= r-1) && (i < r) && (i > 0)
            xi = [ xi; [j i]];
        elseif (0 <= (i+j)) && ((i+j) <= r)
            xb = [ xb; [j i]];
        end
    end
end

figure(1)
clf
hold on
scatter(xb(:, 1), xb(:, 2), 450, '.')
scatter(xi(:, 1), xi(:, 2), 450, '.')
ylabel('i (delay)')
xlabel('j (time range index)')
title('Time Delay Shifting Indices', 'Fontsize', 14)
ylim([0, 6.1])
hold off
% % M = ones(length(i), length(r));
% scr = i' + j;
% 
% % M(i'+j > r) = 0;
% % M(i'+j < -1) = 0;
% M = (scr <= r) .* (scr >= 0);
% spy(M);
% 
