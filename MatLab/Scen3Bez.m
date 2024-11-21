clc;
clear;
close all;
load Scen2Data.mat;
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
AC2_lon = AC2(1, :);
AC2_lat = AC2(2, :);
AC2_hdg = AC2(4, :);
AC3_lon = AC3(1, :);
AC3_lat = AC3(2, :);
AC3_hdg = AC3(4, :);
AC4_lon = AC4(1, :);
AC4_lat = AC4(2, :);
AC4_hdg = AC4(4, :);
EM0_lon = EM0(1, :);
EM0_lat = EM0(2, :);

AC0_lon2 = AC0(1, 1:1:550);
AC0_lat2 = AC0(2, 1:1:550);
AC1_lon2 = AC1(1, 1:1:605);
AC1_lat2 = AC1(2, 1:1:605);

box('on')
ax = gca;
ax.FontWeight = 'bold';
ax.FontName = 'Times';
set(gcf, 'Position', [700 350 1000 1000])

axis('equal')
hold on
r = 0.0004144027532220207;
grid()

plot([EM0_lon(1), EM0_lon(700)], [EM0_lat(1) EM0_lat(700)], 'LineWidth', 3, 'Color', 'red')
% plot([-82.1994, -82.197], [39.599988478353495, 39.6025], 'LineStyle', '--', 'Color', 'k')
rectangle('Position', [EM0_lon(i)-r, EM0_lat(i)-r, 2*r, 2*r], 'Curvature', [1, 1], 'EdgeColor','r', 'LineWidth',2)
scatter(-82.20158972364715, 39.599988478353495, 1000, 'MarkerEdgeColor', on)
text(-82.20425, 39.60175, sprintf('Fleet\nGoal'), FontWeight='bold', FontSize=18, FontName='Times')
scatter(-82.1994, 39.599988478353495, 1000, 'MarkerEdgeColor', 'red')
text(-82.2, 39.6035, sprintf('Emergency\nVehicle Goal'), FontWeight='bold', FontSize=18, FontName='Times')
plot([EM0_lon(352), -82.197], [39.599988478353495, 39.6025], LineStyle='--', Color = 'black', LineWidth=2)


plot(AC0_lon2, AC0_lat2, 'LineWidth', 3, 'Color', 'blue')
plot(AC1_lon2, AC1_lat2, 'LineWidth', 2, 'Color', 'green')

plot([EM0_lon(1), EM0_lon(1)], [EM0_lat(333), EM0_lat(400)], 'LineWidth', 3, 'Color', 'red')


quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(AC1_lon(i), AC1_lat(i), cosd(90-AC1_hdg(i)), sind(90-AC1_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(EM0_lon(i), EM0_lat(i), 1, 0, 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)

scatter(AC0_lon(i), AC0_lat(i), 100, 'o', 'filled', 'blue')
scatter(EM0_lon(i), EM0_lat(i), 100, '^',  'filled', 'red')
scatter(AC1_lon(i), AC1_lat(i), 100, 'o', 'filled', 'green')

text(AC0_lon(i)+0.00055, AC0_lat(i), sprintf('AX0, t_1'), 'FontSize', 18, FontName='Times', FontWeight= 'bold')
text(AC1_lon(i)+0.00055, AC1_lat(i), sprintf('AX1, t_1'), 'FontSize', 18, FontName='Times', FontWeight= 'bold')
text(EM0_lon(i)-0.00075, EM0_lat(i)-.0025, sprintf('EX0, t_1\nGoal TOA\n34.2s'), 'FontSize', 18, FontName='Times', FontWeight= 'bold')



i =243;
rectangle('Position', [EM0_lon(i)-r, EM0_lat(i)-r, 2*r, 2*r], 'Curvature', [1, 1], 'EdgeColor','r', 'LineWidth',2)


quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(AC1_lon(i), AC1_lat(i), cosd(90-AC1_hdg(i)), sind(90-AC1_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(EM0_lon(i), EM0_lat(i), 1, 0, 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)

scatter(AC0_lon(i), AC0_lat(i), 100, 'o', 'filled', 'blue')
scatter(EM0_lon(i), EM0_lat(i), 100, '^',  'filled', 'red')
scatter(AC1_lon(i), AC1_lat(i), 100, 'o', 'filled', 'green')

text(AC0_lon(i)+0.00045, AC0_lat(i)-0.00035, sprintf('AX0, t_2'), 'FontSize', 18, FontName='Times', FontWeight= 'bold')
text(AC1_lon(i)+0.00075, AC1_lat(i), sprintf('AX1, t_2'), 'FontSize', 18, FontName='Times', FontWeight= 'bold')
text(EM0_lon(i)-0.004, EM0_lat(i)-.0025, sprintf('EX0, t_2\nGoal TOA\n10s'), 'FontSize', 18, FontName='Times', FontWeight= 'bold')



% 
i = 299;
rectangle('Position', [EM0_lon(i)-r, EM0_lat(i)-r, 2*r, 2*r], 'Curvature', [1, 1], 'EdgeColor','r', 'LineWidth',2)


quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(AC1_lon(i), AC1_lat(i), cosd(90-AC1_hdg(i)), sind(90-AC1_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(EM0_lon(i), EM0_lat(i), 1, 0, 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)

scatter(AC0_lon(i), AC0_lat(i), 100, 'o', 'filled', 'blue')
scatter(EM0_lon(i), EM0_lat(i), 100, '^',  'filled', 'red')
scatter(AC1_lon(i), AC1_lat(i), 100, 'o', 'filled', 'green')


% plot([AC0_lon(i), AC1_lon(i)], [AC0_lat(i), AC1_lat(i)], 'LineStyle', '--', Color='black', LineWidth=2)
% plot([-82.202, -82.201], [39.5885, 39.595], 'LineStyle', '--', 'Color', 'k')
% 

% text((-82.2+ -82.2021)/2 - 0.0005, (EM0_lat(i)-0.0005), '46m', 'FontSize', 18, FontName='Times', FontWeight= 'bold')
% plot([-82.2, -82.2021], [AC4_lat(i), ACe4_lat(i)], 'LineStyle', '--', 'Color', 'black', 'LineWidth', 3)
% plot([-82.2, -82.2021], [EM0_lat(i), EM0_lat(i)], 'LineStyle', '--', 'Color', 'black', 'LineWidth', 3)
% text(-82.205, 39.5875, sprintf('Sequenc\nViolation\nat t_3'), 'FontSize', 18, FontName='Times', FontWeight= 'bold')

text(AC0_lon(i)+1.5*r, AC0_lat(i), sprintf('\nAX0, t_3'), 'FontSize', 18, FontName='Times', FontWeight= 'bold')
text(AC1_lon(i)+0.00045, AC1_lat(i)+0.0003, sprintf('AX1, t_3'), 'FontSize', 18, FontName='Times', FontWeight= 'bold')
text(EM0_lon(i)-0.003, EM0_lat(i)-.0025, sprintf('EX0, t_3\nGoal TOA\n4.4s'), 'FontSize', 18, FontName='Times', FontWeight= 'bold')

i = 352;
rectangle('Position', [EM0_lon(i)-r, EM0_lat(i)-r, 2*r, 2*r], 'Curvature', [1, 1], 'EdgeColor','r', 'LineWidth',2)


quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(AC1_lon(i), AC1_lat(i), cosd(90-AC1_hdg(i)), sind(90-AC1_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
scatter(AC0_lon(i), AC0_lat(i), 100, 'o', 'filled', 'blue')
scatter(EM0_lon(i), EM0_lat(i), 100, '^',  'filled', 'red')
scatter(AC1_lon(i), AC1_lat(i), 100, 'o', 'filled', 'green')


text(AC0_lon(i)+1*r, AC0_lat(i), sprintf('AX0, t_4'), 'FontSize', 18, FontName='Times', FontWeight= 'bold')
text(AC1_lon(i)+1.5*r, AC1_lat(i)+1*r, sprintf('AX1, t_4'), 'FontSize', 18, FontName='Times', FontWeight= 'bold')
text(EM0_lon(i)-5*r, EM0_lat(i)-3*r, sprintf('\nEX0\nt_4, t_5'), 'FontSize', 15, FontName='Times', FontWeight= 'bold')

% 
% 
i = 450;
rectangle('Position', [EM0_lon(i)-r, EM0_lat(i)-r, 2*r, 2*r], 'Curvature', [1, 1], 'EdgeColor','r', 'LineWidth',3)

quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(AC1_lon(i), AC1_lat(i), cosd(90-AC1_hdg(i)), sind(90-AC1_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)


scatter(AC0_lon(i), AC0_lat(i), 100, 'o', 'filled', 'blue')
scatter(EM0_lon(i), EM0_lat(i), 100, '^',  'filled', 'red')
scatter(AC1_lon(i), AC1_lat(i), 100, 'o', 'filled', 'green')
% plot([AC0_lon(i), AC1_lon(i)], [AC0_lat(i), AC1_lat(i)], 'LineStyle', '--', Color='black', LineWidth=2)


% text((-82.2+ -82.2021)/2 - 0.0005, (EM0_lat(i)-0.0005), '46m', 'FontSize', 18, FontName='Times', FontWeight= 'bold')
% plot([-82.2, -82.2021], [AC4_lat(i), AC4_lat(i)], 'LineStyle', '--', 'Color', 'black', 'LineWidth', 3)
% plot([-82.2, -82.2021], [EM0_lat(i), EM0_lat(i)], 'LineStyle', '--', 'Color', 'black', 'LineWidth', 3)

plot([AC0_lon(i)-2*r, AC0_lon(i)], [AC0_lat(i)+7*r, AC0_lat(i)], 'LineStyle', '--', 'Color', 'k')
plot([AC1_lon(i)-5*r, AC1_lon(i)], [AC1_lat(i)+6*r, AC1_lat(i)], 'LineStyle', '--', 'Color', 'k')

text(AC0_lon(i)-8.5*r, AC0_lat(i)+7*r, sprintf('AX0, t_5'), 'FontSize', 18, FontName='Times', FontWeight= 'bold')
text(AC1_lon(i)-9*r, AC1_lat(i)+7*r, sprintf('AX1, t_5'), 'FontSize', 18, FontName='Times', FontWeight= 'bold')

%
xlim([EM0_lon(1)-0.002, EM0_lon(701)+0.0055])
ylim([AC1_lat(1)-0.002, EM0_lat(400)+0.005])
ylabel('Latitude', FontName='Times')
xlabel('Longitude', FontName='Times')


