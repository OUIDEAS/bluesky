clc;
clear;
close all;
load Scen1DataHold.mat;
brown = [0.75, 0.16, 0.16];
orange = [0.4940 0.1840 0.5560];
purple = [0.900 0.50 0.10];
i = 224
if i == 1
    t = 0
else
    t = i/10
end

set(gcf, 'position', [1000 500 400 600])
box('on')
ax = gca;
ax.FontWeight = 'bold';
ax.FontName = 'Times';
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
% Wayx = Waypts(1, :);
% Wayy = Waypts(2, :);
% Entryx = Entry(1, :);
% Entryy = Entry(2, :);
% Exitx = Exit(1, :);
% Exity = Exit(2, :);

AC0_lon2 = AC0(1, 167:1:329);
AC0_lat2 = AC0(2, 167:1:329);
AC1_lon2 = AC1(1, 226:1:380);
AC1_lat2 = AC1(2, 226:1:380);
AC2_lon2 = AC2(1, 278:1:440);
AC2_lat2 = AC2(2, 278:1:440);
AC3_lon2 = AC3(1, 333:1:495);
AC3_lat2 = AC3(2, 333:1:495);
AC4_lon2 = AC4(1, 387:1:551);
AC4_lat2 = AC4(2, 387:1:551);


% subplot(2, 3, 1)
axis('equal')
hold on
r = 0.0004144027532220207
grid()

plot([-82.2 -82.2], [EM0_lat(i)+0.01 AC4_lat(i)-0.01], 'LineWidth', 3, 'Color', 'black')
rectangle('Position', [EM0_lon(i)-r, EM0_lat(i)-r, 2*r, 2*r], 'Curvature', [1, 1], 'EdgeColor','r', 'LineWidth',2)
plot(AC0_lon2, AC0_lat2, 'LineWidth', 3, 'Color', 'blue')
% plot(AC1_lon2, AC1_lat2, 'LineWidth', 2, 'Color', 'red')
% plot(AC2_lon2, AC2_lat2, 'LineWidth', 2, 'Color', 'green')
% plot(AC3_lon2, AC3_lat2, 'LineWidth', 2, 'Color', orange)
% plot(AC4_lon2, AC4_lat2, 'LineWidth', 2, 'Color', purple)
quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(AC1_lon(i), AC1_lat(i), cosd(90-AC1_hdg(i)), sind(90-AC1_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(AC2_lon(i), AC2_lat(i), cosd(90-AC2_hdg(i)), sind(90-AC2_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(AC3_lon(i), AC3_lat(i), cosd(90-AC3_hdg(i)), sind(90-AC3_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(AC4_lon(i), AC4_lat(i), cosd(90-AC4_hdg(i)), sind(90-AC4_hdg(i)), 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)
quiver(EM0_lon(i), EM0_lat(i), 0, 1, 0.001, 'Color', purple, 'LineWidth', 4, 'MaxHeadSize',4)

scatter(AC0_lon(i), AC0_lat(i), 150, 'o', 'filled', 'blue')
scatter(EM0_lon(i), EM0_lat(i), 150, '^',  'filled', 'red')
scatter(AC1_lon(i), AC1_lat(i), 150, 'o', 'filled', 'green')
scatter(AC2_lon(i), AC2_lat(i), 150, 'o', 'filled', Color=orange)
scatter(AC3_lon(i), AC3_lat(i), 150, 'o', 'filled', Color=purple)
scatter(AC4_lon(i), AC4_lat(i), 150, 'o', 'filled', Color=brown)

% plot(AC1_lon, AC1_lat)
% plot(AC2_lon, AC2_lat)
% plot(AC3_lon, AC3_lat)
% plot(AC4_lon, AC4_lat)

text(-82.204, 39.434, 't = 22.3s', FontWeight='bold', FontSize=18, FontAngle='italic', FontName='Times');
text(-82.20455, 39.4275, sprintf('AC0\nHolding\nPattern'), FontWeight='bold', FontSize=14, FontName='Times')
% grid()
% legend('Aircraft 0', 'Emergency Vehicle','Aircraft 1', 'Aircraft 2', 'Aircraft 3', 'Aircraft 4', 'FontSize', 20)
% text(-8.2, 39.46, 't = 0')
% plot(Wayx, Wayy, 'magenta', 'LineWidth', 2)
% plot(Entryx, Entryy, 'cyan', 'LineWidth', 2)
% plot(Exitx, Exity, 'cyan', 'LineWidth', 2)
plot([AC0_lon(i), AC1_lon(i)], [AC0_lat(i), AC1_lat(i)], 'LineStyle', '--', Color='black', LineWidth=2)

% legend('Aircraft 0', 'Emergency Vehicle','Aircraft 1', 'Aircraft 2', 'Aircraft 3', 'Aircraft 4', 'Partial Bezier Curve', 'FontSize', 20)
% 
% text(EM0_lon(i)+r/2, EM0_lat(i)-r, 'Emergency Vehicle Safety Radius', FontSize=8)
% annotation('textbox', [0.26, 0.7, 0.2, 0.1], 'String', ['After 13.57 seconds Aircraft 0 is within 10 seconds of the emergency vehicle, so a Bezier path' ...
%     'is generated ahead of Aircraft 0'], ...
%     'HorizontalAlignment','center', VerticalAlignment='middle')
% % 
% text(-82.2+r, (EM0_lat(i)+AC0_lat(i))/2+0.000015, '46m Separation', 'FontSize', 10)
% plot([-82.2, -82.1985], [AC0_lat(i), AC0_lat(i)], 'LineStyle', '--', Color='black')
% plot([-82.2, -82.1985], [EM0_lat(i), EM0_lat(i)], 'LineStyle', '--', Color='black')

% text(-82.1985, 39.432, sprintf('Entry Path'),  FontWeight='bold', FontSize=12)
% text(-82.1975, 39.4345, sprintf('Bezier Path'),  FontWeight='bold', FontSize=12)
% text(-82.1985, 39.4365, sprintf('Exit Path'),  FontWeight='bold', FontSize=12)

text(AC0_lon(i)-0.001955, AC0_lat(i)+0.001, sprintf('AC0 at t'), 'FontSize', 14, FontName='Times', FontWeight= 'bold')
text(AC1_lon(i)+0.00055, AC1_lat(i), sprintf('AC1 at t\nViolates\nSequence'), 'FontSize', 14, FontName='Times', FontWeight= 'bold')
text(AC2_lon(i)+0.00055, AC2_lat(i), sprintf('AC2 at t'), 'FontSize', 14, FontName='Times', FontWeight= 'bold')
text(AC3_lon(i)+0.00055, AC3_lat(i), sprintf('AC3 at t'), 'FontSize', 14, FontName='Times', FontWeight= 'bold')
text(AC4_lon(i)+0.00055, AC4_lat(i), sprintf('AC4 at t'), 'FontSize', 14, FontName='Times', FontWeight= 'bold')
text(EM0_lon(i)+r, EM0_lat(i)-.0001, sprintf('EM0 at t'), 'FontSize', 14, FontName='Times', FontWeight= 'bold')
% text(-82.2-10*r, EM0_lat(i)+0.0005, sprintf('\nEmergency Vehicle\n46m Safety Radius'), 'FontSize', 11, FontName='Times', FontWeight= 'bold')
ylim([39.4204-0.0015, EM0_lat(i)+0.0032])
ylabel('Latitude', FontName='Times')
xlabel('Longitude', FontName='Times')
% ZoomPlot2()
% [p,z] = zoomPlot(AC0_lon, AC0_lat, [-82.2 -82.2], [.6 0.65 .25 .25], [1 4])
% hold on \nViolates Sequence And\nInitiates Holding Pattern