% Data
clc;
clear;
close all;
set(gcf, 'Position', [1000 500 750 700])
x_labels = {'AX0', 'AX1'};%, 'AC2', 'AC3', 'AC4'};
ETA = [
    %Scen1
    % 348.5412648256665, 350.4037558669337, 368.8912648256665; 
    % 353.9630462925526, 353.9630462925526, 374.55304629255255;]; 
    %Scen2
    % 348.5412648256665, 350.40354299908756, 368.89126482566644;
    % 353.9630462925526, 356.00908149753207, 374.3130462925526;];
    %Scen3
    34.13963219284771, 43.234879048968935, 54.48963219284771;
    39.55005658715413, 47.10494388093571, 60.140056587154135;];

% Extract columns for plotting
eta_0 = ETA(:, 1);
eta_1 = ETA(:, 2);
eta_2 = ETA(:, 3);

% Bar width
bar_width = 0.35;

% X positions for the bars
x = 1:length(x_labels);

% Create a figure
% figure;5


% Plot the bars
hold on;
% bar0 = bar(x - bar_width, eta_0, bar_width, 'FaceColor', 'green'); % Green
bar1 = bar(x-0.175, (eta_2-eta_0), bar_width, 'FaceColor', 'red'); % Blue
bar2 = bar(x-0.175 + bar_width, eta_1-eta_0, bar_width, 'FaceColor', 'blue');   % Red 1.3275
plot([0, x(1)-bar_width/2], [eta_2(1)-eta_0(1), eta_2(1)-eta_0(1)], 'linestyle', '--', 'color', 'k', 'linewidth', 2)
plot([0, x(2)-bar_width/2], [eta_2(2)-eta_0(2), eta_2(2)-eta_0(2)], 'linestyle', '--', 'color', 'k', 'linewidth', 2)
plot([0, x(1)+bar_width/2], [eta_1(1)-eta_0(1), eta_1(1)-eta_0(1)], 'linestyle', '--', 'color', 'k', 'linewidth', 2)
% text(x(2)+bar_width/1.25, 2.5, sprintf('No\ndelay\nabosrbed'), 'FontSize', 25, FontName='Times', FontWeight= 'bold', HorizontalAlignment= 'center', VerticalAlignment = 'middle')
plot([0, x(2)+bar_width/2], [eta_1(2)-eta_0(2), eta_1(2)-eta_0(2)], 'linestyle', '--', 'color', 'k', 'linewidth', 2)
plot([0, 0], [eta_1(2)-eta_0(2), eta_1(2)-eta_0(2)], 'linestyle', '--', 'color', 'k', 'linewidth', 2)


% for i = 1:length(x_labels)
%     % Calculate percent changes
%     pct_change_1 = ((eta_1(i) - eta_0(i)) / eta_0(i)) * 100;
%     pct_change_2 = ((eta_2(i) - eta_0(i)) / eta_0(i)) * 100;
%     del_1 = eta_1(i) - eta_0(i);
%     del_2 = eta_2(i) - eta_0(i);
% 
%     % Display percent change between eta_2 and eta_0
%     text(x(i) + 1.125*bar_width,  del_1+1, sprintf('%.2fs%', del_1), ...
%         'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 18); %eta_1(i)
%     % if i == 1
%     %     text(x(i)-.45 + bar_width/3-0.1, eta_1(i) + 2.5, sprintf('%% Change: %.2f%%', pct_change_1), ...
%     %         'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 18);
%     % else
%     %     text(x(i) + bar_width/3-0.15, eta_1(i) + 2.5, sprintf('%.2f%%', pct_change_1), ...
%     %         'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 18);
%     % end
% 
%     % Display percent change between eta_2 and eta_0
%     % text(x(i) + bar_width, eta_2(i) + 2.5, sprintf('%.2f%%', pct_change_2), ...
%     %     'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 18);
%     if i == 1
%         text(x(i) + bar_width/3-0.1, del_2 +1, sprintf('%.2fs% s', del_2), ...
%             'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 18);%eta_2(i)
%     else
%         text(x(i) + bar_width/3-0.05, del_2 +1, sprintf('%.2fs% s', del_2), ...
%             'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 18);
%     end
% end

% Set x-ticks and labels
set(gca, 'XTick', x, 'XTickLabel', x_labels);
grid()
ax = gca;
ax.FontSize = 18;
% ylim([min(eta_0) - 10, max(eta_2) + 10]);
ylim([0, 25]);
xlim([x(1)-2*bar_width, x(2)+2*bar_width])
% Labeling
xlabel('Aircraft', 'FontSize', 18);
ylabel('Aircraft Delay (s)', 'FontSize', 18);
title('Scenario 3 Aircraft Delays', 'FontSize', 18);
legend({'Holding Pattern Alternate Maneuver Delay', 'Bezier Alternate Maneuever Delay'}, 'Location','northwest', 'FontSize', 13);%'Nominal ETA', 


% Display the plot
hold off;
% Data

%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Percent Change Graph  %
%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% clc;
% clear;
% close all;
% set(gcf, 'Position', [1000 500 750 700])
% x_labels = {'AX0', 'AX1'};%, 'AC2', 'AC3', 'AC4'};
% ETA = [
%     348.5412648256665, 350.4037558669337, 368.8912648256665; 
%     353.9630462925526, 353.9630462925526, 374.55304629255255;]; 
%     % 44.96320583097026, 52.518093123953214, 65.79320583097027; 
%     % 50.37820184019237, 57.933089135728956, 71.44820184019237; 
%     % 55.79450715868693, 63.349394453424885, 77.10450715868694];
% 
% % Extract columns for plotting
% eta_0 = ETA(:, 1);
% eta_1 = ETA(:, 2);
% eta_2 = ETA(:, 3);
% 
% % Bar width
% bar_width = 0.35;
% 
% % X positions for the bars
% x = 1:length(x_labels);
% 
% % Create a figure
% % figure;
% 
% 
% % Plot the bars
% hold on;
% % bar0 = bar(x - bar_width, eta_0, bar_width, 'FaceColor', 'green'); % Green
% bar1 = bar(x(1), 100*((eta_2(1)-eta_0(1))/eta_0(1)), bar_width, 'FaceColor', 'red'); % Blue
% bar2 = bar(x(1) + bar_width, 100*((eta_1(1)-eta_0(1))/eta_0(1)), bar_width, 'FaceColor', 'blue');   % Red 1.3275
% bar3 = bar(x(2), 100*((eta_2(2)-eta_0(2))/eta_0(2)), bar_width, 'FaceColor', 'red'); % Blue
% bar4 = bar(x(2) + bar_width, 100*((eta_1(2)-eta_0(2))/eta_0(2)), bar_width, 'FaceColor', 'blue');   % Red 1.3275
% hold off;
% for i = 1:length(x_labels)
% %     % Calculate percent changes
%     pct_change_1 = ((eta_1(i) - eta_0(i)) / eta_0(i)) * 100;
%     pct_change_2 = ((eta_2(i) - eta_0(i)) / eta_0(i)) * 100;
%     del_1 = eta_1(i) - eta_0(i);
%     del_2 = eta_2(i) - eta_0(i);
% % 
% %     % Display percent change between eta_2 and eta_0
% %     text(x(i) + 1.125*bar_width,  del_1+1, sprintf('%.2fs%', del_1), ...
% %         'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 18); %eta_1(i)
% %     % if i == 1
%         text(x(i) + bar_width+0.02, pct_change_1+.2, sprintf('%.2f%%', pct_change_1), ...
%             'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 18);
% %     % else
% %     %     text(x(i) + bar_width/3-0.15, eta_1(i) + 2.5, sprintf('%.2f%%', pct_change_1), ...
% %     %         'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 18);
% %     % end
% % 
% %     % Display percent change between eta_2 and eta_0
%     text(x(i), pct_change_2 + .2, sprintf('%.2f%%', pct_change_2), ...
%         'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 18);
% %     if i == 1
% %         text(x(i) + bar_width/3-0.1, del_2 +1, sprintf('%.2fs% s', del_2), ...
% %             'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 18);%eta_2(i)
% %     else
% %         text(x(i) + bar_width/3-0.05, del_2 +1, sprintf('%.2fs% s', del_2), ...
% %             'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 18);
% %     end
% end
% 
% % Set x-ticks and labels
% set(gca, 'XTick', x, 'XTickLabel', x_labels);
% grid()
% ax = gca;
% ax.FontSize = 18;
% % ylim([min(eta_0) - 10, max(eta_2) + 10]);
% ylim([0, 7]);
% % Labeling
% xlabel('Aircraft', 'FontSize', 18);
% ylabel('% Change in ETA', 'FontSize', 18);
% title('Scenario 1 Aircraft % Change in ETA', 'FontSize', 18);
% legend({'Holding Pattern Alternate Maneuver', 'Bezier Alternate Maneuever'}, 'Location','northwest', 'FontSize', 13);%'Nominal ETA', 
% 
% 
% % Display the plot
% % hold off;
