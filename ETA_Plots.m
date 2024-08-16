% Data
clc;
clear;
close all;
set(gcf, 'Position', [1000 500 750 700])
x_labels = {'AC0', 'AC1'};%, 'AC2', 'AC3', 'AC4'};
ETA = [
    34.13963219284771, 43.234879048968935, 54.48963219284771; 
    39.55005658715413, 47.10494388093571, 60.140056587154135;]; 
    % 44.96320583097026, 52.518093123953214, 65.79320583097027; 
    % 50.37820184019237, 57.933089135728956, 71.44820184019237; 
    % 55.79450715868693, 63.349394453424885, 77.10450715868694];

% Extract columns for plotting
eta_0 = ETA(:, 1);
eta_1 = ETA(:, 2);
eta_2 = ETA(:, 3);

% Bar width
bar_width = 0.35;

% X positions for the bars
x = 1:length(x_labels);

% Create a figure
% figure;


% Plot the bars
hold on;
bar0 = bar(x - bar_width/3, eta_0, bar_width, 'FaceColor', 'green'); % Green
bar1 = bar(x + bar_width/3, eta_1, bar_width, 'FaceColor', 'blue'); % Blue
bar2 = bar(x + 1.1*bar_width, eta_2, bar_width, 'FaceColor', 'red');   % Red

for i = 1:length(x_labels)
    % Calculate percent changes
    pct_change_1 = ((eta_1(i) - eta_0(i)) / eta_0(i)) * 100;
    pct_change_2 = ((eta_2(i) - eta_0(i)) / eta_0(i)) * 100;
    del_1 = eta_1(i) - eta_0(i);
    del_2 = eta_2(i) - eta_0(i);
    
    % Display percent change between eta_2 and eta_0
    text(x(i) + bar_width, eta_2(i) +1, sprintf('%.2fs%', del_2), ...
        'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold');
    if i == 1
        text(x(i)-.25 + bar_width/3-0.1, eta_1(i) + 2.5, sprintf('Percent Change: %.2f%%', pct_change_1), ...
            'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold');
    else
        text(x(i) + bar_width/3-0.05, eta_1(i) + 2.5, sprintf('%.2f%%', pct_change_1), ...
            'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold');
    end

    % Display percent change between eta_2 and eta_0
    text(x(i) + bar_width, eta_2(i) + 2.5, sprintf('%.2f%%', pct_change_2), ...
        'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold');
    if i == 1
        text(x(i)-.19 + bar_width/3-0.1, eta_1(i) +1, sprintf('Total Delay: %.2fs% s', del_1), ...
            'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold');
    else
        text(x(i) + bar_width/3-0.05, eta_1(i) +1, sprintf('%.2fs% s', del_1), ...
            'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold');
    end
end

% Set x-ticks and labels
set(gca, 'XTick', x, 'XTickLabel', x_labels);
grid()
ylim([min(eta_0) - 10, max(eta_2) + 10]);
% Labeling
xlabel('Aircraft');
ylabel('ETA Values (s)');
title('Scenario 3 Aircraft ETA');
legend({'Nominal ETA','Bezier Alternate Maneuever ETA', 'Holding Pattern Alternate Maneuver ETA'}, 'Location','northeast');


% Display the plot
hold off;
