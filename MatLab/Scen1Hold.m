clc;
clear;
close all;
load Scen1DataHold.mat;
brown = [0.75, 0.16, 0.16];
orange = [0.4940 0.1840 0.5560];
purple = [0.900 0.50 0.10];
c1 = [0.635 0.078 0.184];
c2 = [0 0.447 0.741];
on = [.5 0 0.5];
i = 1;

AC0_lon = AC0(1, :);
AC0_lat = AC0(2, :);
AC0_hdg = AC0(4, :);
AC1_lon = AC1(1, :);
AC1_lat = AC1(2, :);
AC1_hdg = AC1(4, :);
EM0_lon = EM0(1, :);
EM0_lat = EM0(2, :);

AC0_lon2 = AC0(1, 1:1:561);
AC0_lat2 = AC0(2, 1:1:561);
AC1_lon2 = AC1(1, 1:1:561);
AC1_lat2 = AC1(2, 1:1:561);


% subplot(3, 3, 1);

box('on')
ax = gca;
ax.FontWeight = 'bold';
ax.FontName = 'Times';
set(gcf, 'Position', [700 350 1550 700])
% 
% subplot(1, 4, 1);
t = tiledlayout(1,4, 'TileSpacing', 'compact', 'Padding', 'compact'); % 1x4 grid with minimal spacing

ax = nexttile;
ax.FontWeight = 'bold';
% text(-82.2075, 39.4475, 'a_0 Path', FontWeight='bold', FontSize=18)
axis('equal')
hold on
r = 0.0004144027532220207;
grid()
i=1;
plot(AC0_lon2, AC0_lat2, 'LineWidth', 3, 'Color', 'blue')
quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
scatter(AC0_lon(i), AC0_lat(i), 100, 'o', 'filled', 'blue')
text(AC0_lon(i)+1.5*r, AC0_lat(i), sprintf('a_0, t_0'), 'FontSize', 18, FontWeight= 'bold')

i =167;
quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
scatter(AC0_lon(i), AC0_lat(i), 100, 'o', 'filled', 'blue')
text(AC0_lon(i)+1.5*r, AC0_lat(i), sprintf('t_1'), 'FontSize', 18, FontWeight= 'bold')

% 
i = 224;
quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
scatter(AC0_lon(i), AC0_lat(i), 100, 'o', 'filled', 'blue')
text(AC0_lon(i)-6*r, AC0_lat(i), sprintf('t_2'), 'FontSize', 18, FontWeight= 'bold')

% plot([AC1_lon(i), AC1_lon(i)+3.5*r], [AC1_lat(i), AC1_lat(i)+3.75*r], 'LineStyle','-', 'Color','k')
% plot([AC0_lon(i), AC1_lon(i)], [AC0_lat(i), AC1_lat(i)], 'LineWidth', 2, 'LineStyle', '--', 'Color', 'k')
% text(-82.209, 39.435, sprintf('Sequence\nViolation, t_2'), 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times')
% plot([-82.206, -82.20125], [39.434, 39.4285], 'LineStyle', '-', 'Color', 'k')
% plot([EM0_lon(i), EM0_lon(i)+3.5*r], [EM0_lat(i), EM0_lat(i)+2*r], 'LineStyle','-', 'Color','k')
% 
i = 330;
quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
scatter(AC0_lon(i), AC0_lat(i), 100, 'o', 'filled', 'blue')
text(AC0_lon(i)+1.5*r, AC0_lat(i), sprintf('t_3'), 'FontSize', 18, FontWeight= 'bold')

% plot([AC0_lon(i), AC0_lon(i)+1.5*r], [AC0_lat(i), AC0_lat(i)-4.5*r], 'Color', 'k')
% plot([AC1_lon(i), AC0_lon(i)], [AC1_lat(i), AC0_lat(i)], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'k')
% text(-82.209, 39.4225, sprintf('Sequence\nRestored, t4'), 'FontSize',18, 'FontWeight','bold', 'FontName','Times')
% plot([-82.2025, -82.20125], [39.422, 39.426], 'Color', 'k')

% 
% 
i = 530;
quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
scatter(AC0_lon(i), AC0_lat(i), 100, 'o', 'filled', 'blue')
text(AC0_lon(i)+1.5*r, AC0_lat(i), sprintf('t_4'), 'FontSize', 18, FontWeight= 'bold')

% plot([AC1_lon(i), AC1_lon(i)+3.5*r], [AC1_lat(i), AC1_lat(i)+2*r], 'Color', 'k')

xlim([EM0_lon(1)-0.01, EM0_lon(560)+0.01])
ylim([AC1_lat(1)-0.0015, EM0_lat(560)+0.0001])
% ylabel('Latitude', FontName='Times')
% xlabel('Longitude', FontName='Times')
xlabel('a_0 Path', FontWeight='bold', FontSize=18)

% % hold on \nViolates Sequence And\nInitiates Holding Pattern
% % \nViolates\nSequence \nReturns\nTo Path\nSequence Is\nResumed
% 
hold off
nexttile;
axis('equal')
hold on
r = 0.0004144027532220207;
grid()
i=1;
% text(-82.2075, 39.4475, 'a_1 Path', FontWeight='bold', FontSize=18)

plot(AC1_lon2, AC1_lat2, 'LineWidth', 2, 'Color', 'green')
quiver(AC1_lon(i), AC1_lat(i), cosd(90-AC1_hdg(i)), sind(90-AC1_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
scatter(AC1_lon(i), AC1_lat(i), 100, 'o', 'filled', 'green')
text(AC1_lon(i)+1.5*r, AC1_lat(i), sprintf('a_1, t_0'), 'FontSize', 18, FontWeight= 'bold')

i =167;
quiver(AC1_lon(i), AC1_lat(i), cosd(90-AC1_hdg(i)), sind(90-AC1_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
scatter(AC1_lon(i), AC1_lat(i), 100, 'o', 'filled', 'green')
text(AC1_lon(i)+1.5*r, AC1_lat(i), sprintf('t_1'), 'FontSize', 18, FontWeight= 'bold')

% plot([EM0_lon(i), EM0_lon(i)+1.5*r], [EM0_lat(i), EM0_lat(i)-0.0005], 'LineStyle', '-', 'Color', 'k')


% 
i = 224;
quiver(AC1_lon(i), AC1_lat(i), cosd(90-AC1_hdg(i)), sind(90-AC1_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
scatter(AC1_lon(i), AC1_lat(i), 100, 'o', 'filled', 'green')
text(AC1_lon(i)+1.5*r, AC1_lat(i), sprintf('t_2'), 'FontSize', 18, FontWeight= 'bold')

% plot([AC1_lon(i), AC1_lon(i)+3.5*r], [AC1_lat(i), AC1_lat(i)+3.75*r], 'LineStyle','-', 'Color','k')
% plot([AC0_lon(i), AC1_lon(i)], [AC0_lat(i), AC1_lat(i)], 'LineWidth', 2, 'LineStyle', '--', 'Color', 'k')
% text(-82.209, 39.435, sprintf('Sequence\nViolation, t_2'), 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times')
% plot([-82.206, -82.20125], [39.434, 39.4285], 'LineStyle', '-', 'Color', 'k')
% plot([EM0_lon(i), EM0_lon(i)+3.5*r], [EM0_lat(i), EM0_lat(i)+2*r], 'LineStyle','-', 'Color','k')
% 
i = 330;
quiver(AC1_lon(i), AC1_lat(i), cosd(90-AC1_hdg(i)), sind(90-AC1_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
scatter(AC1_lon(i), AC1_lat(i), 100, 'o', 'filled', 'green')
text(AC1_lon(i)-6*r, AC1_lat(i), sprintf('t_3'), 'FontSize', 18, FontWeight= 'bold')
% 
% plot([AC0_lon(i), AC0_lon(i)+1.5*r], [AC0_lat(i), AC0_lat(i)-4.5*r], 'Color', 'k')
% plot([AC1_lon(i), AC0_lon(i)], [AC1_lat(i), AC0_lat(i)], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'k')
% text(-82.209, 39.4225, sprintf('Sequence\nRestored, t4'), 'FontSize',18, 'FontWeight','bold', 'FontName','Times')
% plot([-82.2025, -82.20125], [39.422, 39.426], 'Color', 'k')

% 
% 
i = 530;
quiver(AC1_lon(i), AC1_lat(i), cosd(90-AC1_hdg(i)), sind(90-AC1_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
scatter(AC1_lon(i), AC1_lat(i), 100, 'o', 'filled', 'green')
text(AC1_lon(i)+1.5*r, AC1_lat(i), sprintf('t_4'), 'FontSize', 18, FontWeight= 'bold')

% plot([AC1_lon(i), AC1_lon(i)+3.5*r], [AC1_lat(i), AC1_lat(i)+2*r], 'Color', 'k')

xlim([EM0_lon(1)-0.01, EM0_lon(560)+0.01])
ylim([AC1_lat(1)-0.0015, EM0_lat(560)+0.0001])
% ylabel('Latitude', FontName='Times')
% xlabel('Longitude', FontName='Times')
set(gca, 'XTickLabel', [], 'YTickLabel', []);
xlabel('a_1 Path', FontWeight='bold', FontSize=18)

hold off
% % hold on \nViolates Sequence And\nInitiates Holding Pattern
% % \nViolates\nSequence \nReturns\nTo Path\nSequence Is\nResumed
% 

nexttile;
axis('equal')
hold on
r = 0.0004144027532220207;
grid()
i=1;
% text(-82.2075, 39.4475, 'e_0 Path', FontWeight='bold', FontSize=18)

plot([EM0_lon(1), EM0_lon(400)], [EM0_lat(1) EM0_lat(561)], 'LineWidth', 3, 'Color', 'red')
plot([EM0_lon(1), EM0_lon(1)], [EM0_lat(333), EM0_lat(400)], 'LineWidth', 3, 'Color', 'red')

quiver(EM0_lon(i), EM0_lat(i), 0, 1, 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
scatter(EM0_lon(i), EM0_lat(i), 100, '^',  'filled', 'red')
text(EM0_lon(i)+1.5*r, EM0_lat(i), sprintf('e_0, t_0'), 'FontSize', 18, FontWeight= 'bold')

i =167;
rectangle('Position', [EM0_lon(i)-r, EM0_lat(i)-r, 2*r, 2*r], 'Curvature', [1, 1], 'EdgeColor','r', 'LineWidth',2)

quiver(EM0_lon(i), EM0_lat(i), 0, 1, 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
scatter(EM0_lon(i), EM0_lat(i), 100, '^',  'filled', 'red')
text(EM0_lon(i)+1.5*r, EM0_lat(i), sprintf('t_1'), 'FontSize', 18, FontWeight= 'bold')

% plot([EM0_lon(i), EM0_lon(i)+1.5*r], [EM0_lat(i), EM0_lat(i)-0.0005], 'LineStyle', '-', 'Color', 'k')


% 
i = 224;
rectangle('Position', [EM0_lon(i)-r, EM0_lat(i)-r, 2*r, 2*r], 'Curvature', [1, 1], 'EdgeColor','r', 'LineWidth',2)

quiver(EM0_lon(i), EM0_lat(i), 0, 1, 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
scatter(EM0_lon(i), EM0_lat(i), 100, '^',  'filled', 'red')
text(EM0_lon(i)+1.5*r, EM0_lat(i), sprintf('t_2'), 'FontSize', 18, FontWeight= 'bold')

% plot([AC1_lon(i), AC1_lon(i)+3.5*r], [AC1_lat(i), AC1_lat(i)+3.75*r], 'LineStyle','-', 'Color','k')
% plot([AC0_lon(i), AC1_lon(i)], [AC0_lat(i), AC1_lat(i)], 'LineWidth', 2, 'LineStyle', '--', 'Color', 'k')
% text(-82.209, 39.435, sprintf('Sequence\nViolation, t_2'), 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times')
% plot([-82.206, -82.20125], [39.434, 39.4285], 'LineStyle', '-', 'Color', 'k')
% plot([EM0_lon(i), EM0_lon(i)+3.5*r], [EM0_lat(i), EM0_lat(i)+2*r], 'LineStyle','-', 'Color','k')
% % 
i = 330;
rectangle('Position', [EM0_lon(i)-r, EM0_lat(i)-r, 2*r, 2*r], 'Curvature', [1, 1], 'EdgeColor','r', 'LineWidth',2)
quiver(EM0_lon(i), EM0_lat(i), 0, 1, 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
scatter(EM0_lon(i), EM0_lat(i), 100, '^',  'filled', 'red')
text(EM0_lon(i)+1.5*r, EM0_lat(i), sprintf('t_3'), 'FontSize', 18, FontWeight= 'bold')
% 
% plot([AC0_lon(i), AC0_lon(i)+1.5*r], [AC0_lat(i), AC0_lat(i)-4.5*r], 'Color', 'k')
% plot([AC1_lon(i), AC0_lon(i)], [AC1_lat(i), AC0_lat(i)], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'k')
% text(-82.209, 39.4225, sprintf('Sequence\nRestored, t4'), 'FontSize',18, 'FontWeight','bold', 'FontName','Times')
% plot([-82.2025, -82.20125], [39.422, 39.426], 'Color', 'k')

% 
% 
i = 530;
rectangle('Position', [EM0_lon(i)-r, EM0_lat(i)-r, 2*r, 2*r], 'Curvature', [1, 1], 'EdgeColor','r', 'LineWidth',3)
quiver(EM0_lon(i), EM0_lat(i), 0, 1, 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
scatter(EM0_lon(i), EM0_lat(i), 100, '^',  'filled', 'red')
text(EM0_lon(i)+1.5*r, EM0_lat(i), sprintf('t_4'), 'FontSize', 18, FontWeight='bold')

% plot([AC1_lon(i), AC1_lon(i)+3.5*r], [AC1_lat(i), AC1_lat(i)+2*r], 'Color', 'k')

xlim([EM0_lon(1)-0.01, EM0_lon(560)+0.01])
ylim([AC1_lat(1)-0.0015, EM0_lat(560)+0.0001])
% ylabel('Latitude', FontName='Times')
% xlabel('Longitude', FontName='Times')
set(gca, 'XTickLabel', [], 'YTickLabel', []);
xlabel('e_0 Path', FontWeight='bold', FontSize=18)

% % hold on \nViolates Sequence And\nInitiates Holding Pattern
% % \nViolates\nSequence \nReturns\nTo Path\nSequence Is\nResumed
% 


hold off

nexttile;
axis('equal')
hold on
% text(-82.2085, 39.4475, 'Total Path', FontWeight='bold', FontSize=18)

r = 0.0004144027532220207;
grid()
i=1;
plot([EM0_lon(1), EM0_lon(400)], [EM0_lat(1) EM0_lat(561)], 'LineWidth', 3, 'Color', 'red')

plot(AC0_lon2, AC0_lat2, 'LineWidth', 3, 'Color', 'blue')
plot(AC1_lon2, AC1_lat2, 'LineWidth', 2, 'Color', 'green')
plot([EM0_lon(1), EM0_lon(1)], [EM0_lat(333), EM0_lat(400)], 'LineWidth', 3, 'Color', 'red')



quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(AC1_lon(i), AC1_lat(i), cosd(90-AC1_hdg(i)), sind(90-AC1_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(EM0_lon(i), EM0_lat(i), 0, 1, 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)

scatter(AC0_lon(i), AC0_lat(i), 100, 'o', 'filled', 'blue')
scatter(EM0_lon(i), EM0_lat(i), 100, '^',  'filled', 'red')
scatter(AC1_lon(i), AC1_lat(i), 100, 'o', 'filled', 'green')

text(AC0_lon(i)+1.5*r, AC0_lat(i), sprintf('a_0, t_0'), 'FontSize', 18, FontWeight= 'bold')
text(AC1_lon(i)+1.5*r, AC1_lat(i), sprintf('a_1, t_0'), 'FontSize', 18, FontWeight= 'bold')
text(EM0_lon(i)-11*r, EM0_lat(i), sprintf('e_0, t_0'), 'FontSize', 18, FontWeight= 'bold')




i =167;
rectangle('Position', [EM0_lon(i)-r, EM0_lat(i)-r, 2*r, 2*r], 'Curvature', [1, 1], 'EdgeColor','r', 'LineWidth',2)

quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(AC1_lon(i), AC1_lat(i), cosd(90-AC1_hdg(i)), sind(90-AC1_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(EM0_lon(i), EM0_lat(i), 0, 1, 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)

scatter(AC0_lon(i), AC0_lat(i), 100, 'o', 'filled', 'blue')
scatter(EM0_lon(i), EM0_lat(i), 100, '^',  'filled', 'red')
scatter(AC1_lon(i), AC1_lat(i), 100, 'o', 'filled', 'green')

text(AC0_lon(i)+1.5*r, AC0_lat(i)+0.00015, sprintf('a_0, t_1'), 'FontSize', 18, FontWeight= 'bold')
text(AC1_lon(i)+1.5*r, AC1_lat(i)-0.001, sprintf('a_1, t_1'), 'FontSize', 18, FontWeight= 'bold')
text(EM0_lon(i)+1.5*r, EM0_lat(i)-0.0005, sprintf('\ne_0, t_1'), 'FontSize', 18, FontWeight= 'bold')

plot([EM0_lon(i), EM0_lon(i)+1.5*r], [EM0_lat(i), EM0_lat(i)-0.0005], 'LineStyle', '-', 'Color', 'k')


% 
i = 224;
rectangle('Position', [EM0_lon(i)-r, EM0_lat(i)-r, 2*r, 2*r], 'Curvature', [1, 1], 'EdgeColor','r', 'LineWidth',2)

quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(AC1_lon(i), AC1_lat(i), cosd(90-AC1_hdg(i)), sind(90-AC1_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(EM0_lon(i), EM0_lat(i), 0, 1, 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)

scatter(AC0_lon(i), AC0_lat(i), 100, 'o', 'filled', 'blue')
scatter(EM0_lon(i), EM0_lat(i), 100, '^',  'filled', 'red')
scatter(AC1_lon(i), AC1_lat(i), 100, 'o', 'filled', 'green')

text(AC0_lon(i)-11*r, AC0_lat(i), sprintf('a_0, t_2'), 'FontSize', 18, FontWeight= 'bold')
text(AC1_lon(i)+1.5*r, AC1_lat(i)+4*r, sprintf('a_1, t_2'), 'FontSize', 18, FontWeight= 'bold')
text(EM0_lon(i)+1.5*r, EM0_lat(i)+2*r, sprintf('e_0, t_2'), 'FontSize', 18, FontWeight= 'bold')

plot([AC1_lon(i), AC1_lon(i)+3.5*r], [AC1_lat(i), AC1_lat(i)+3.75*r], 'LineStyle','-', 'Color','k')
plot([AC0_lon(i), AC1_lon(i)], [AC0_lat(i), AC1_lat(i)], 'LineWidth', 2, 'LineStyle', '--', 'Color', 'k')
text(-82.209, 39.435, sprintf('Sequence\nViolation, t_2'), 'FontSize', 16, 'FontWeight', 'bold')
plot([-82.206, -82.20125], [39.434, 39.4285], 'LineStyle', '-', 'Color', 'k')
plot([EM0_lon(i), EM0_lon(i)+3.5*r], [EM0_lat(i), EM0_lat(i)+2*r], 'LineStyle','-', 'Color','k')
% 
i = 330;
rectangle('Position', [EM0_lon(i)-r, EM0_lat(i)-r, 2*r, 2*r], 'Curvature', [1, 1], 'EdgeColor','r', 'LineWidth',2)


quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(AC1_lon(i), AC1_lat(i), cosd(90-AC1_hdg(i)), sind(90-AC1_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(EM0_lon(i), EM0_lat(i), 0, 1, 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)

scatter(AC0_lon(i), AC0_lat(i), 100, 'o', 'filled', 'blue')
scatter(EM0_lon(i), EM0_lat(i), 100, '^',  'filled', 'red')
scatter(AC1_lon(i), AC1_lat(i), 100, 'o', 'filled', 'green')


text(AC0_lon(i)+1.5*r, AC0_lat(i)-6*r, sprintf('a_0, t_3'), 'FontSize', 18, FontWeight= 'bold')
text(AC1_lon(i)-11*r, AC1_lat(i), sprintf('a_1, t_3'), 'FontSize', 18, FontWeight= 'bold')
text(EM0_lon(i)+1.5*r, EM0_lat(i)+2*r, sprintf('e_0, t_3'), 'FontSize', 18, FontWeight= 'bold')

plot([AC0_lon(i), AC0_lon(i)+1.5*r], [AC0_lat(i), AC0_lat(i)-4.5*r], 'Color', 'k')
plot([AC1_lon(i), AC0_lon(i)], [AC1_lat(i), AC0_lat(i)], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'k')
text(-82.209, 39.4225, sprintf('Sequence\nRestored, t4'), 'FontSize',16, 'FontWeight','bold')
plot([-82.2025, -82.20125], [39.422, 39.426], 'Color', 'k')

% 
% 
i = 530;
rectangle('Position', [EM0_lon(i)-r, EM0_lat(i)-r, 2*r, 2*r], 'Curvature', [1, 1], 'EdgeColor','r', 'LineWidth',3)

quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(AC1_lon(i), AC1_lat(i), cosd(90-AC1_hdg(i)), sind(90-AC1_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(EM0_lon(i), EM0_lat(i), 0, 1, 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)

scatter(AC0_lon(i), AC0_lat(i), 100, 'o', 'filled', 'blue')
scatter(EM0_lon(i), EM0_lat(i), 100, '^',  'filled', 'red')
scatter(AC1_lon(i), AC1_lat(i), 100, 'o', 'filled', 'green')


text(AC0_lon(i)+1.5*r, AC0_lat(i), sprintf('a_0, t_4'), 'FontSize', 18, FontWeight= 'bold')
text(AC1_lon(i)+1.5*r, AC1_lat(i)+2*r, sprintf('a_1, t_4'), 'FontSize', 18, FontWeight= 'bold')
text(EM0_lon(i)+1.5*r, EM0_lat(i), sprintf('e_0, t_4'), 'FontSize', 18, FontWeight='bold')

plot([AC1_lon(i), AC1_lon(i)+3.5*r], [AC1_lat(i), AC1_lat(i)+2*r], 'Color', 'k')

xlim([EM0_lon(1)-0.01, EM0_lon(560)+0.01])
ylim([AC1_lat(1)-0.0015, EM0_lat(560)+0.0001])
% ylabel('Latitude', FontName='Times')
% xlabel('Longitude', FontName='Times')
set(gca, 'XTickLabel', [], 'YTickLabel', []);
xlabel('Combined Path', FontWeight='bold', FontSize=18)

hold off
% % hold on \nViolates Sequence And\nInitiates Holding Pattern
% % \nViolates\nSequence \nReturns\nTo Path\nSequence Is\nResumed
% 
xlabel(t, 'Longitude', 'FontSize', 18, FontWeight = 'bold');
ylabel(t, 'Latitude', 'FontSize', 18, FontWeight = 'bold');