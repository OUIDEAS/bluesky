clc;
clear;
close all;
box('on');
set(gcf, 'Position', [1000 500 800 750]);
brown = [200 100 16]/256;
orange = [0.8500 0.3250 0.0980];
purple = [0.4940 0.1840 0.5560];
forest_green = [34, 139, 34] / 255;


% load P5XBez.mat;
hold on;
% plot([Bez1X, Bez2X], [Bez1Y, Bez2Y], 'Color', 'magenta', 'LineWidth', 2);
% plot(, Bez2Y, 'Color', 'magenta', 'LineWidth', 2);
% 
load 1XBez.mat;

plot( [Bez1Y, Bez2Y],-[Bez1X, Bez2X], 'Color', 'blue', 'LineWidth', 2);
% plot(Bez2X, Bez2Y, 'Color', 'blue', 'LineWidth', 2);
% 
load 1P5XBez.mat;

plot( [Bez1Y, Bez2Y],-[Bez1X, Bez2X], 'Color', 'green', 'LineWidth', 2);
% plot(Bez2X, Bez2Y, 'Color', 'green', 'LineWidth', 2);

load 2XBez.mat;

plot( [Bez1Y, Bez2Y],-[Bez1X, Bez2X], 'Color', purple, 'LineWidth', 2);

load 1XHold.mat
plot( [top_y, straight_y, bot_y, straight_y],-[turn_x, right_x, turn_x2, left_x, ], LineWidth=2, Color='blue', LineStyle='--')
% plot(, , LineWidth=2, Color='black')
% plot(, , LineWidth=2, Color='black')
% plot(, , LineWidth=2, Color='black')

load 1P5XHold.mat
plot([top_y, straight_y, bot_y, straight_y],-[turn_x, right_x, turn_x2, left_x, ],  LineWidth=2, Color='green', LineStyle='--')

load 2XHold.mat
plot( [top_y, straight_y, bot_y, straight_y],-[turn_x, right_x, turn_x2, left_x, ], LineWidth=2, Color=purple, LineStyle='--')
% plot(Bez2X, Bez2Y, 'Color', purple, 'LineWidth', 2);
% 
% load 3XBez.mat;
% 
% plot( [Bez1Y, Bez2Y],-[Bez1X, Bez2X], 'Color', 'magenta', 'LineWidth', 2);
% % plot(Bez2X, Bez2Y, 'Color', 'black', 'LineWidth', 2);
% 
% load 4XBez.mat;
% 
% plot( [Bez1Y, Bez2Y],-[Bez1X, Bez2X], 'Color', 'cyan', 'LineWidth', 2);
% % plot(Bez2X, Bez2Y, 'Color', 'cyan', 'LineWidth', 2);


load Corridors.mat
plot( [B1_NodeY(1), B1_NodeY(2), B2_NodeY(3)],-[corridor, corridor, corridor], 'LineStyle', '--', 'color', orange, LineWidth = 3);
% plot([op_corridor-228, op_corrridor-228, op_corridor-228], [B1_NodeY(1), B1_NodeY(2), B2_NodeY(3)], 'LineStyle', '--', 'color', 'black');
plot( [koz_top, koz_bot,koz_bot, koz_bot, koz_top, koz_top],-[koz_x, koz_x,B1_NodeX(1), koz_x, koz_x, B1_NodeX(1)], color = 'red', linestyle = '--', LineWidth = 3);
% plot([], [], color = 'red', linestyle = '--', LineWidth = 3);
% plot([], [], color = 'red', linestyle = '--', LineWidth = 3);
axis('equal');
grid();
ylim([-300, 200]);
xlim([0, 1100]);
% title('Bezier Curves with Different Requested TOAs', FontSize=14);
ylabel('Y(m)', FontSize=14);
xlabel('X(m)', FontSize=14);
% show();

% load P5XHold.mat
% plot([turn_x, right_x, turn_x2, left_x, ], [top_y, straight_y, bot_y, straight_y], LineWidth=2, Color='black')
% plot([turn_x, right_x], [top_y, straight_y], LineWidth=2, Color='black')


% load 3XHold.mat
% plot( [top_y, straight_y, bot_y, straight_y],-[turn_x, right_x, turn_x2, left_x, ], LineWidth=2, Color='magenta', LineStyle='--')
% 
% load 4XHold.mat
% plot( [top_y, straight_y, bot_y, straight_y],-[turn_x, right_x, turn_x2, left_x, ], LineWidth=2, Color='cyan', LineStyle='--')

plot([-600, 0, 1100], [0, 0, 0], LineWidth=2, color = forest_green, HandleVisibility= 'on')

scatter(B2_NodeY(3), B2_NodeX(3), 100, MarkerFaceColor='black', MarkerEdgeColor='black')

scatter(B1_NodeY(1), B1_NodeX(1), 100, MarkerFaceColor='green', MarkerEdgeColor='green')

legend({'5.3s', '7.95s', '10.6s', '5.3s', '7.95s', '10.6s', 'Flight Corridor Bound', 'Emergency Vehicle' + string(newline)  + 'Clearance Area',  'Nominal Path', 'Goal Point', 'Start Point'}, 'Location','north west', 'FontSize', 13, 'NumColumns', 4);

plot( [B1_NodeY(3), B1_NodeY(3), B1_NodeY(3)],-[0, 150, 250], LineStyle='--', color = 'black', LineWidth=2, HandleVisibility= 'off')
text( 475,-100, sprintf('Curve 1'), FontSize=10)
text( 850,-150, sprintf('Curve 2'), FontSize=10)
text(200, -250, sprintf('Holding Patterns'), FontSize=10)

%'15.9s', '21.2s' '15.9s', '21.2s',