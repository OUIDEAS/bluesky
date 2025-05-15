clc;
clear;
close all;
load crash.mat;
brown = [0.75, 0.16, 0.16];
orange = [0.4940 0.1840 0.5560];
purple = [0.900 0.50 0.10];
c1 = [0.635 0.078 0.184];
c2 = [0 0.447 0.741];
on = [.5 0 0.5];

AC0_lon = AC0(1, :);
AC0_lat = AC0(2, :);
AC0_hdg = AC0(4, :);

EM0_lon = EM0(1, :);
EM0_lat = EM0(2, :);

AC0_lon2 = AC0(1, 1:1:701);
AC0_lat2 = AC0(2, 1:1:701);

box('on')
ax = gca;
ax.FontWeight = 'bold';
ax.FontName = 'Times';
set(gcf, 'Position', [700 350 1250 700])
% 
% subplot(1, 4, 1);
t = tiledlayout(1,3, 'TileSpacing', 'compact', 'Padding', 'compact'); % 1x4 grid with minimal spacing
ax = nexttile;
ax.FontWeight = 'bold';
axis('equal')
hold on
r = 0.0004144027532220207/5;
grid()
% text(-82.2045, 39.443, 'a_0 Path', 'FontSize', 18, FontWeight= 'bold')
i = 1;

% plot(AC0_lon2, AC0_lat2, 'LineWidth', 3, 'Color', 'blue')
plot([AC0_lon(1), AC0_lon(207)], [AC0_lat(1) AC0_lat(207)], 'LineWidth', 3, 'Color', 'blue')
quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
scatter(AC0_lon(i), AC0_lat(i), 100, 'o', 'filled', 'blue')
text(AC0_lon(i)+5*r, AC0_lat(i), sprintf('a_0, t_0'), 'FontSize', 18, FontWeight= 'bold')

i =100;
quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)

scatter(AC0_lon(i), AC0_lat(i), 100, 'o', 'filled', 'blue')

text(AC0_lon(i)+1.5*5*r, AC0_lat(i), sprintf('t_1'), 'FontSize', 18, FontWeight= 'bold')

i  = 208;
quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color' , purple, 'LineWidth', 4, 'MaxHeadSize',4)

scatter(AC0_lon(i), AC0_lat(i), 100, 'o', 'filled', 'blue')
text(AC0_lon(i)+1.5*5*r, AC0_lat(i), sprintf('t_2'), 'FontSize', 18, FontWeight= 'bold')

% i  = 400;
% quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
% 
% scatter(AC0_lon(i), AC0_lat(i), 100, 'o', 'filled', 'blue')
% 
% text(AC0_lon(i)+1.5*5*r, AC0_lat(i), sprintf('t_3'), 'FontSize', 18, FontWeight= 'bold')
% i  = 400;
% 
% quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
% scatter(AC0_lon(i), AC0_lat(i), 100, 'o', 'filled', 'blue')
% text(AC0_lon(i)+1.5*5*r, AC0_lat(i), sprintf('t_4'), 'FontSize', 18, FontWeight= 'bold')
% xlim([EM0_lon(1)-0.0045, EM0_lon(207)+0.0045])
% ylim([AC1_lat(1)-0.0015, EM0_lat(401)+0.0001])
ylim([39.5875, 39.6025])% ylabel('Latitude', FontName='Times')
% xlabel('Longitude', FontName='Times')
xlabel('a_0 Path', FontWeight='bold', FontSize=18)

hold off

nexttile;
% 
axis('equal')
hold on
r = 0.0004144027532220207/5;
grid()

i=1;
plot([EM0_lon(1), EM0_lon(207)], [EM0_lat(1) EM0_lat(207)], 'LineWidth', 3, 'Color', 'red')
% plot([EM0_lon(1), EM0_lon(1)], [EM0_lat(333), EM0_lat(207)], 'LineWidth', 3, 'Color', 'red')
quiver(EM0_lon(i), EM0_lat(i), 1, 0, 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
scatter(EM0_lon(i), EM0_lat(i), 100, '^',  'filled', 'red')
text(EM0_lon(i)-7.5*5*r, EM0_lat(i)-5*5*r, sprintf('e_0,t_0'), 'FontSize', 18, FontWeight= 'bold')

i =100;
rectangle('Position', [EM0_lon(i)-r, EM0_lat(i)-r, 2*r, 2*r], 'Curvature', [1, 1], 'EdgeColor','r', 'LineWidth',2)
quiver(EM0_lon(i), EM0_lat(i), 1, 0, 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
scatter(EM0_lon(i), EM0_lat(i), 100, '^',  'filled', 'red')
text(EM0_lon(i), EM0_lat(i) - 5*5*r, sprintf('t_1'), 'FontSize', 18, FontWeight= 'bold')

i  = 208;
rectangle('Position', [EM0_lon(i)-r, EM0_lat(i)-r, 2*r, 2*r], 'Curvature', [1, 1], 'EdgeColor','r', 'LineWidth',2)
quiver(EM0_lon(i), EM0_lat(i), 1, 0, 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
scatter(EM0_lon(i), EM0_lat(i), 100, '^',  'filled', 'red')
text(EM0_lon(i), EM0_lat(i) - 5*5*r, sprintf('t_2'), 'FontSize', 18, FontWeight= 'bold')

% i  = 400;
% rectangle('Position', [EM0_lon(i)-r, EM0_lat(i)-r, 2*5*r, 2*5*r], 'Curvature', [1, 1], 'EdgeColor','r', 'LineWidth',2)
% quiver(EM0_lon(i), EM0_lat(i), 1, 0, 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
% scatter(EM0_lon(i), EM0_lat(i), 100, '^',  'filled', 'red')
% text(EM0_lon(i), EM0_lat(i) - 5*5*r, sprintf('t_3'), 'FontSize', 18, FontWeight= 'bold')

% i  = 400;
% rectangle('Position', [EM0_lon(i)-r, EM0_lat(i)-r, 2*5*r, 2*5*r], 'Curvature', [1, 1], 'EdgeColor','r', 'LineWidth',3)
% quiver(EM0_lon(i), EM0_lat(i), 1, 0, 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
% scatter(EM0_lon(i), EM0_lat(i), 100, '^',  'filled', 'red')
% text(EM0_lon(i), EM0_lat(i) - 5*5*r, sprintf('t_4'), 'FontSize', 18, FontWeight='bold')

% set(gca, 'XTickLabel', [], 'YTickLabel', []);
xlim([EM0_lon(1)-0.005, EM0_lon(207)+0.005])
% % ylim([AC1_lat(1)-0.0015, EM0_lat(401)+0.0001])
ylim([39.5875, 39.605])
xlabel('e_0 Path', FontWeight='bold', FontSize=18)

% ylabel('Latitude', FontName='Times')
% xlabel('Longitude', FontName='Times')
hold off

nexttile;
axis('equal')
hold on
r = 0.0004144027532220207/5;
grid on;
% text(-82.2045, 39.4435, 'Total Path', 'FontSize', 18, FontWeight= 'bold')

i=1;
plot([EM0_lon(1), EM0_lon(207)], [EM0_lat(1) EM0_lat(207)], 'LineWidth', 3, 'Color', 'red')
plot([AC0_lon(1), AC0_lon(207)], [AC0_lat(1) AC0_lat(207)], 'LineWidth', 3, 'Color', 'blue')
% plot(AC0_lon2, AC0_lat2, 'LineWidth', 3, 'Color', 'blue')
% plot(AC1_lon2, AC1_lat2, 'LineWidth', 2, 'Color', 'green')
plot([EM0_lon(1), EM0_lon(1)], [EM0_lat(333), EM0_lat(207)], 'LineWidth', 3, 'Color', 'red')

quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
% quiver(AC1_lon(i), AC1_lat(i), cosd(90-AC1_hdg(i)), sind(90-AC1_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(EM0_lon(i), EM0_lat(i), 1, 0, 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)

scatter(AC0_lon(i), AC0_lat(i), 100, 'o', 'filled', 'blue')
scatter(EM0_lon(i), EM0_lat(i), 100, '^',  'filled', 'red')
% scatter(AC1_lon(i), AC1_lat(i), 100, 'o', 'filled', 'green')

text(AC0_lon(i)+5*r, AC0_lat(i), sprintf('a_0, t_0'), 'FontSize', 17.5, FontWeight= 'bold')
% text(AC1_lon(i)+3*5*r, AC1_lat(i), sprintf('a_1, t_0'), 'FontSize', 18, FontWeight= 'bold')
text(EM0_lon(i)-7.5*5*r, EM0_lat(i)-5*5*r, sprintf('e_0,t_0'), 'FontSize', 18, FontWeight= 'bold')

i =100;
rectangle('Position', [EM0_lon(i)-r, EM0_lat(i)-r, 2*r, 2*r], 'Curvature', [1, 1], 'EdgeColor','r', 'LineWidth',2)

quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
% quiver(AC1_lon(i), AC1_lat(i), cosd(90-AC1_hdg(i)), sind(90-AC1_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(EM0_lon(i), EM0_lat(i), 1, 0, 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)

scatter(AC0_lon(i), AC0_lat(i), 100, 'o', 'filled', 'blue')
scatter(EM0_lon(i), EM0_lat(i), 100, '^',  'filled', 'red')
% scatter(AC1_lon(i), AC1_lat(i), 100, 'o', 'filled', 'green')

text(AC0_lon(i)+1.5*5*r, AC0_lat(i), sprintf('t_1'), 'FontSize', 18, FontWeight= 'bold')
% text(AC1_lon(i)+3*5*r, AC1_lat(i), sprintf('a_1, t_1'), 'FontSize', 18, FontWeight= 'bold')
text(EM0_lon(i), EM0_lat(i)-5*5*r, sprintf('t_1'), 'FontSize', 18, FontWeight= 'bold')

% plot([EM0_lon(i), EM0_lon(i)+3*5*r], [EM0_lat(i), EM0_lat(i)-0.0005], 'LineStyle', '-', 'Color', 'k')


i  = 208;
rectangle('Position', [EM0_lon(i)-r, EM0_lat(i)-r, 2*r, 2*r], 'Curvature', [1, 1], 'EdgeColor','r', 'LineWidth',2)

quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
% quiver(AC1_lon(i), AC1_lat(i), cosd(90-AC1_hdg(i)), sind(90-AC1_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(EM0_lon(i), EM0_lat(i), 1, 0, 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)

scatter(AC0_lon(i), AC0_lat(i), 100, 'o', 'filled', 'blue')
scatter(EM0_lon(i), EM0_lat(i), 100, '^',  'filled', 'red')
% scatter(AC1_lon(i), AC1_lat(i), 100, 'o', 'filled', 'green')

text(AC0_lon(i)+1.5*5*r, AC0_lat(i)+7*r, sprintf('t_2'), 'FontSize', 18, FontWeight= 'bold')
% text(AC1_lon(i)+1.5*5*r, AC1_lat(i), sprintf('a_1, t_2'), 'FontSize', 18, FontWeight= 'bold')
% text(EM0_lon(i), EM0_lat(i)-5*5*r, sprintf('t_2'), 'FontSize', 18, FontWeight= 'bold')

% plot([AC1_lon(i), AC1_lon(i)+1.5*5*r], [AC1_lat(i), AC1_lat(i)+3.75*5*r], 'LineStyle','-', 'Color','k')
% plot([AC0_lon(i), AC1_lon(i)], [AC0_lat(i), AC1_lat(i)], 'LineWidth', 2, 'LineStyle', '--', 'Color', 'k')
% text(-82.209, 39.435, sprintf('Sequence\nViolation, t_2'), 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times')
% plot([-82.206, -82.20125], [39.434, 39.4285], 'LineStyle', '-', 'Color', 'k')
% plot([EM0_lon(i), EM0_lon(i)+1.5*5*r], [EM0_lat(i), EM0_lat(i)+2*5*r], 'LineStyle','-', 'Color','k')
% 
% i  = 400;
% rectangle('Position', [EM0_lon(i)-r, EM0_lat(i)-r, 2*5*r, 2*5*r], 'Curvature', [1, 1], 'EdgeColor','r', 'LineWidth',2)
% 
% 
% quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
% % quiver(AC1_lon(i), AC1_lat(i), cosd(90-AC1_hdg(i)), sind(90-AC1_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
% quiver(EM0_lon(i), EM0_lat(i), 1, 0, 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
% 
% scatter(AC0_lon(i), AC0_lat(i), 100, 'o', 'filled', 'blue')
% scatter(EM0_lon(i), EM0_lat(i), 100, '^',  'filled', 'red')
% % scatter(AC1_lon(i), AC1_lat(i), 100, 'o', 'filled', 'green')
% 
% 
% text(AC0_lon(i)+1.5*5*r, AC0_lat(i), sprintf('t_3'), 'FontSize', 18, FontWeight= 'bold')
% % text(AC1_lon(i)-7.5*5*r, AC1_lat(i), sprintf('a_1, t_3'), 'FontSize', 18, FontWeight= 'bold')
% text(EM0_lon(i), EM0_lat(i)-5*5*r, sprintf('t_3'), 'FontSize', 18, FontWeight= 'bold')

% plot([AC0_lon(i), AC0_lon(i)+3*5*r], [AC0_lat(i), AC0_lat(i)-4.5*5*r], 'Color', 'k')
% plot([AC1_lon(i), AC0_lon(i)], [AC1_lat(i), AC0_lat(i)], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'k')
% text(-82.209, 39.4225, sprintf('Sequence\nRestored, t4'), 'FontSize',18, 'FontWeight','bold', 'FontName','Times')
% plot([-82.2025, -82.20125], [39.422, 39.426], 'Color', 'k')

% 
% 
% i  = 400;
% rectangle('Position', [EM0_lon(i)-r, EM0_lat(i)-r, 2*5*r, 2*5*r], 'Curvature', [1, 1], 'EdgeColor','r', 'LineWidth',3)
% 
% quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
% % quiver(AC1_lon(i), AC1_lat(i), cosd(90-AC1_hdg(i)), sind(90-AC1_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
% quiver(EM0_lon(i), EM0_lat(i), 1, 0, 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
% 
% scatter(AC0_lon(i), AC0_lat(i), 100, 'o', 'filled', 'blue')
% scatter(EM0_lon(i), EM0_lat(i), 100, '^',  'filled', 'red')
% % scatter(AC1_lon(i), AC1_lat(i), 100, 'o', 'filled', 'green')
% % 
% 
% text(AC0_lon(i)-7.5*5*r, AC0_lat(i), sprintf('a_0, t_4'), 'FontSize', 18, FontWeight= 'bold')
% % text(AC1_lon(i)+1.5*5*r, AC1_lat(i), sprintf('a_1, t_4'), 'FontSize', 18, FontWeight= 'bold')
% text(EM0_lon(i), EM0_lat(i) - 5*5*r, sprintf('e_0,t_4'), 'FontSize', 18, FontWeight='bold')

% plot([AC1_lon(i), AC1_lon(i)+1.5*5*r], [AC1_lat(i), AC1_lat(i)+2*5*r], 'Color', 'k')
% ZoomPlot2()
% [p,z] = zoomPlot(AC0_lon, AC0_lat, [-82.2 -82.2], [.6 0.65 .25 .25], [1 4])
xlim([EM0_lon(1)-0.005, EM0_lon(207)+0.005])
% ylim([AC1_lat(1)-0.0015, EM0_lat(401)+0.0001])
ylim([39.5875, 39.605])
% ylabel('Latitude', FontName='Times')
% xlabel('Longitude', FontName='Times')
set(gca, 'XTickLabel', [], 'YTickLabel', []);
xlabel('Combined Path', FontWeight='bold', FontSize=18)

% grid()
hold off
% % hold on \nViolates Sequence And\nInitiates Holding Pattern
% % \nViolates\nSequence \nReturns\nTo Path\nSequence Is\nResumed
% 
xlabel(t, 'Longitude', 'FontSize', 18, FontWeight = 'bold');
ylabel(t, 'Latitude', 'FontSize', 18, FontWeight = 'bold');