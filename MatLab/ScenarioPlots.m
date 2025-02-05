clc;
clear;
close all;
load FullBezR.mat;
brown = [0.75, 0.16, 0.16];
orange = [0.4940 0.1840 0.5560];
purple = [0.900 0.50 0.10];
c1 = [0.635 0.078 0.184];
c2 = [0 0.447 0.741];
on = [.5 0 0.5];
i = 250;

AC0_lon = AC0(1, :);
AC0_lat = AC0(2, :);
AC0_hdg = AC0(4, :);
AC1_lon = AC1(1, :);
AC1_lat = AC1(2, :);
AC1_hdg = AC1(4, :);
% AC2_lon = AC2(1, :);
% AC2_lat = AC2(2, :);
% AC2_hdg = AC2(4, :);
% AC3_lon = AC3(1, :);
% AC3_lat = AC3(2, :);
% AC3_hdg = AC3(4, :);
% AC4_lon = AC4(1, :);
% AC4_lat = AC4(2, :);
% AC4_hdg = AC4(4, :);
EM0_lon = EM0(1, :);
EM0_lat = EM0(2, :);

% AC0_lon2 = AC0(1, 242:1:450);
% AC0_lat2 = AC0(2, 242:1:450);
% AC1_lon2 = AC1(1, 297:1:505);
% AC1_lat2 = AC1(2, 297:1:505);
% AC2_lon2 = AC2(1, 353:1:558);
% AC2_lat2 = AC2(2, 353:1:558);
% AC3_lon2 = AC3(1, 408:1:616);
% AC3_lat2 = AC3(2, 408:1:616);
% AC4_lon2 = AC4(1, 463:1:671);
% AC4_lat2 = AC4(2, 463:1:671);

% figure('Position', [1000 100 1200 1200]);

% subplot(3, 3, 1);

box('on')
ax = gca;
ax.FontWeight = 'bold';
ax.FontName = 'Times';
set(gcf, 'Position', [1000 500 400 600])

axis('equal')
hold on
r = 0.0004144027532220207;
grid()

% plot([-82.2 -82.2], [EM0_lat(i)+0.05 AC4_lat(i)-0.01], 'LineWidth', 3, 'Color', 'black')
rectangle('Position', [EM0_lon(i)-r, EM0_lat(i)-r, 2*r, 2*r], 'Curvature', [1, 1], 'EdgeColor','r', 'LineWidth',2)
% scatter(-82.20158972364715, 39.599988478353495, 500, 'pentagram', 'filled', 'MarkerFaceColor', on)
% text(-82.205, 39.601, sprintf('Fleet\nGoal'), FontWeight='bold', FontSize=11, FontName='Times')
% scatter(-82.1994, 39.599988478353495, 500, 'pentagram', 'filled', 'MarkerFaceColor', 'red')
% text(-82.2, 39.6025, sprintf('Emergency\nVehicle Goal'), FontWeight='bold', FontSize=11, FontName='Times')
% plot([-82.20158972364715,-82.2, -82.2], [39.599988478353495, AC3_lat(408),AC4_lat(1)-0.0032], 'LineWidth', 3,'Color', 'black')
% text(EM0_lon(i)-0.00155, EM0_lat(i)+0.003, sprintf('GOAL ETA:\n0s'), FontWeight='bold', FontSize=14, FontName='Times')
% plot(AC0_lon2, AC0_lat2, 'LineWidth', 3, 'Color', 'blue')
% plot(AC1_lon2, AC1_lat2, 'LineWidth', 2, 'Color', 'green')
% plot(AC2_lon2, AC2_lat2, 'LineWidth', 2, 'Color',  c1)
% plot(AC3_lon2, AC3_lat2, 'LineWidth', 2, 'Color', c2)
% plot(AC4_lon2, AC4_lat2, 'LineWidth', 2, 'Color', brown)

quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(AC1_lon(i), AC1_lat(i), cosd(90-AC1_hdg(i)), sind(90-AC1_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(AC2_lon(i), AC2_lat(i), cosd(90-AC2_hdg(i)), sind(90-AC2_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(AC3_lon(i), AC3_lat(i), cosd(90-AC3_hdg(i)), sind(90-AC3_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(AC4_lon(i), AC4_lat(i), cosd(90-AC4_hdg(i)), sind(90-AC4_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
% quiver(EM0_lon(i), EM0_lat(i), 1, 0, 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)

scatter(AC0_lon(i), AC0_lat(i), 150, 'o', 'filled', 'blue')
scatter(EM0_lon(i), EM0_lat(i), 150, '^',  'filled', 'red')
scatter(AC1_lon(i), AC1_lat(i), 150, 'o', 'filled', 'green')
scatter(AC2_lon(i), AC2_lat(i), 150, 'o', 'filled', 'MarkerEdgeColor', c1, 'MarkerFaceColor', c1)
scatter(AC3_lon(i), AC3_lat(i), 150, 'o', 'filled', 'MarkerEdgeColor', c2, 'MarkerFaceColor', c2)
scatter(AC4_lon(i), AC4_lat(i), 150, 'o', 'filled', 'MarkerEdgeColor', brown, 'MarkerFaceColor', brown)



text(-82.1972, 39.6025, 't = 35.2s', FontWeight='bold', FontSize=18, FontAngle='italic', FontName='Times');

text(-82.20575, 39.593, sprintf('AC2\nHolding\nPattern'), FontWeight='bold', FontSize=10, FontName='Times')
% grid()

plot([AC2_lon(i), AC1_lon(i)], [AC2_lat(i), AC1_lat(i)], 'LineStyle', '--', Color='black', LineWidth=2)
% 

% text((-82.2+ -82.2021)/2 - 0.0005, (EM0_lat(i)-0.0005), '46m', 'FontSize', 14, FontName='Times', FontWeight= 'bold')
% plot([-82.2, -82.2021], [AC4_lat(i), AC4_lat(i)], 'LineStyle', '--', 'Color', 'black', 'LineWidth', 3)
% plot([-82.2, -82.2021], [EM0_lat(i), EM0_lat(i)], 'LineStyle', '--', 'Color', 'black', 'LineWidth', 3)

text(AC0_lon(i)-11.5*r, AC0_lat(i), sprintf('\n\nAC0 at t'), 'FontSize', 14, FontName='Times', FontWeight= 'bold')
text(AC1_lon(i)-11.5*r, AC1_lat(i), sprintf('AC1 at t'), 'FontSize', 14, FontName='Times', FontWeight= 'bold')
text(AC2_lon(i)+0.00055, AC2_lat(i), sprintf('\nAC2 at t\nViolates Sequence'), 'FontSize', 14, FontName='Times', FontWeight= 'bold')
text(AC3_lon(i)+0.00055, AC3_lat(i), sprintf('AC3 at t'), 'FontSize', 14, FontName='Times', FontWeight= 'bold')
text(AC4_lon(i)+0.00055, AC4_lat(i), sprintf('AC4 at t'), 'FontSize', 14, FontName='Times', FontWeight= 'bold')
text(EM0_lon(i)-0.00155, EM0_lat(i)-.0001, sprintf('\n\nEM0 at t'), 'FontSize', 14, FontName='Times', FontWeight= 'bold')
% text(-82.2-12*r, EM0_lat(i)+0.0005, sprintf('\nEmergency Vehicle\n46m Safety Radius'), 'FontSize', 11, FontName='Times', FontWeight= 'bold')
% ylim([39.4204-0.0015, EM0_lat(i)+0.0032])
% ylim([39.58, 39.605])
% xlim([EM0_lon(242)-0.0015, -82.19])
% ZoomPlot2()
ylabel('Latitude', FontName='Times')
xlabel('Longitude', FontName='Times')

% hold on \nViolates Sequence And\nInitiates Holding Pattern
% \nViolates\nSequence \nReturns\nTo Path\nSequence Is\nResumed

