clc; clear; close all;
load tbone_dist.mat;
load crash_dist.mat;
x = 2069;
t = 1:1:x;
% t;
% dist_crash;
hold on;
plot(t/10, dist_crash(1, 1:1:x), linewidth = 3)
x_s = zeros(1,x-1);
for i = 1:1:x-1
    x_s(i) = x/10;
end
y_s = 1:1800/x:1800;
min(dist_evasion)
x = 2200;
t = 1:1:x;
plot(t/10, dist_evasion(1, 1:1:x), linewidth = 3)

% plot(x_s, y_s, linestyle='--', color='black', linewidth = 2)
% plot(206.9, dist_crash(2069), 'x', 'Color', 'red', 'LineWidth', 3, 'MarkerSize', 12)
% text(200, 300, sprintf('Collision\nPoint'), fontweight = 'bold', fontsize = 16)

plot([100, 100], [-50, 1800], linewidth = 2, linestyle = '--', color = 'black')
text(105, 1400, 't_1', fontweight = 'bold', fontsize = 16)
plot([206.8, 206.9], [-50, 1800], linewidth = 2, linestyle = '--', color = 'black')
text(211, 1400, 't_2', fontweight = 'bold', fontsize = 16)
% plot([100, 100], [-50, 1800], linewidth = 2, linestyle = '--', color = 'black')

grid('on')
xt = get(gca, 'XTick');         % Get current x-tick positions
set(gca, 'XTickLabel', xt / 10) % Set the new labels divided by 10
ylabel('Distance Between Aircraft', fontweight = 'bold', FontSize=16)
xlabel('Time (s)', fontweight = 'bold', FontSize=16)
xlim([0, 230])
ylim([-50, 1800])
legend('No Maneuver', 'Bezier Maneuver', fontsize = 10)