clc; clear; close all;

function [curvex,curvey] = quadBezier(P0, P1, P2, points)
    t = linspace(0,1,points);
    curvex = P1(1) + (P0(1)-P1(1))*(1 - t).^2 + (P2(1)-P1(1))*t.^2;
    curvey = P1(2) + (P0(2)-P1(2))*(1 - t).^2 + (P2(2)-P1(2))*t.^2;
end

function [lon_points, lat_points] = Meters_To_WSG84(cx1, cy1, home)
    % Convert position from UTM back to latitude/longitude
    % cx1, cy1: Lists of x and y points in meters
    % home: [latitude, longitude] in WGS84
    % lon_points, lat_points: Separate lists of longitude and latitude points

    % Define the UTM zone 17 projection (WGS84)
    utmproj = projcrs(32617, 'Authority', 'EPSG'); % EPSG 32617: UTM zone 17N

    % Convert home latitude and longitude to UTM
    [homeX, homeY] = projfwd(utmproj, home(1), home(2));

    % Ensure cx1 and cy1 are column vectors
    cx1 = cx1(:);
    cy1 = cy1(:);

    % Initialize output arrays
    n_points = length(cx1);
    lat_points = zeros(n_points, 1);
    lon_points = zeros(n_points, 1);

    % Convert each point back to latitude/longitude
    for i = 1:n_points
        % Add offsets to the home coordinates
        x = homeX + cx1(i);
        y = homeY + cy1(i);

        % Convert from UTM back to lat/lon
        [lat, lon] = projinv(utmproj, x, y);

        % Store results
        lat_points(i) = lat;
        lon_points(i) = lon;
    end
end


path = 'C:/Users/Michael/Desktop/bluesky/MatLab/SQLData/';

data = struct();
xl = zeros([1, 10]);
c = 1;


for sp = 8:8:104
    sp_str = sprintf('sp%.0f', sp); 
    ac_data = load(sprintf('%s%.0f_aircraft.mat', path, sp));
    bez_data = load(sprintf('%s%.0f_bez.mat', path, sp));
    delay_data = load(sprintf('%s%.0f_delay.mat', path, sp));
    dub_data = load(sprintf('%s%.0f_dub.mat', path, sp));
    ev_data = load(sprintf('%s%.0f_ev_specific.mat', path, sp));
    note_data = load(sprintf('%s%.0f_note_events.mat', path, sp));
    state_data = load(sprintf('%s%.0f_state.mat', path, sp));

    data.(sp_str).aircraft = ac_data;
    data.(sp_str).bez = bez_data;
    data.(sp_str).delay = delay_data;
    data.(sp_str).dub = dub_data;
    data.(sp_str).ev_specific = ev_data;
    data.(sp_str).note_events = note_data;
    data.(sp_str).state = state_data;

    xl(c) = sp;
    c = c+1;
end

bar_width = 1;
hold on;

% x = 1:length(xl)
% for i = 1:length(x)
% sp_str = sprintf('sp%.0f', xl(i))
%  delay_value = data.(sp_str).delay.delay_dat.delay;
% 
%     % Plot using x(i) and delay_value
%     plot([0, x(i) + bar_width * i], [delay_value(1), delay_value(1)], 'LineWidth', 2);
% end

num_aircraft = 2; % Assuming 2 aircraft per sp
delay_matrix = zeros(length(xl), num_aircraft);

% Loop through each sp to populate the delay matrix
for i = 1:length(xl)
    sp_str = sprintf('sp%.0f', xl(i)); % Create the dynamic field name

    if data.(sp_str).delay.delay_dat.delay > 0
    delay_values = data.(sp_str).delay.delay_dat.delay; % Extract delay values for this sp
    if xl(i) == 311
        delay_values = delay_values(1:2);
    end
    
    delay_matrix(i, :) = delay_values; 
    end

   % Store the delay values in the matrix
end

% Create bar plot with grouped bars
x = 1:length(xl); % X values corresponding to sp indices
bar(x, delay_matrix, bar_width); % Plot grouped bars

% Customize the plot
xticks(x);
% grid();
xticklabels(arrayfun(@num2str, xl, 'UniformOutput', false)); % Set x-axis labels to sp numbers
xlabel('Spacing Between Fleet Aircraft (m)');
ylabel('Delay (s)');
legend({'AX0', 'AX1'}, 'Location', 'best');
title('Delay for Each Aircraft Across SP Values');
hold off;



% home = [39.49029047106718, -82.2];
% 
% [cx1, cy1] = quadBezier([0, 574.12], [205.6159, 564.0093], [221, 711.62], 100);
% 
% [bezlon1, bezlat1] = Meters_To_WSG84(cx1, cy1, home)
% 
% 
% lon = state_dat.lon(13550:1:17001);
% lat = state_dat.lat(13550:1:17001);
% 
% 
% 
% hold on;
% plot(lon, lat, LineWidth=2);
% plot(bezlon1, bezlat1, LineWidth=2)
% xlim([-82.2-0.003, -82.2+0.003])
% axis('equal')


