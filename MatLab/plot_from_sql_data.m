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


path = 'C:/Users/Michael/Desktop/bluesky/MatLab/SingleHoldSQLData/';

data = struct();
xl = zeros([1, 3]);
mins = zeros([0, 3]);
c = 1;

path2 = 'C:/Users/Michael/Desktop/bluesky/MatLab/SingleBezSQLData/';

data2 = struct();
xl2 = zeros([1, 3]);
mins2 = zeros([0, 3]);
c2 = 1;

path3 = 'C:/Users/Michael/Desktop/bluesky/MatLab/FullHoldSQLData/';

data3 = struct();
xl3 = zeros([1, 3]);
mins3 = zeros([0, 3]);
c3 = 1;
for sp = 24:-8:8
    sp_str = sprintf('sp%.0f', sp); 
    ac_data = load(sprintf('%s%.0f_aircraft.mat', path, sp));
    % bez_data = load(sprintf('%s%.0f_bez.mat', path, sp));
    delay_data = load(sprintf('%s%.0f_delay.mat', path, sp));
    % dub_data = load(sprintf('%s%.0f_dub.mat', path, sp));
    ev_data = load(sprintf('%s%.0f_ev_specific.mat', path, sp));
    note_data = load(sprintf('%s%.0f_note_events.mat', path, sp));
    state_data = load(sprintf('%s%.0f_state.mat', path, sp));
    % dist_data = load(sprintf('%s%.0f_distances.mat', path, sp));

    data.(sp_str).aircraft = ac_data;
    % data.(sp_str).bez = bez_data;
    data.(sp_str).delay = delay_data;
    % data.(sp_str).dub = dub_data;
    data.(sp_str).ev_specific = ev_data;
    data.(sp_str).note_events = note_data;
    data.(sp_str).state = state_data;
    % data.(sp_str).distances = dist_data;

    xl(c) = sp;
    c = c+1;

    ac_data2 = load(sprintf('%s%.0f_aircraft.mat', path2, sp));
    % bez_data2 = load(sprintf('%s%.0f_bez.mat', path2, sp));
    delay_data2 = load(sprintf('%s%.0f_delay.mat', path2, sp));
    % dub_data2 = load(sprintf('%s%.0f_dub.mat', path2, sp));
    ev_data2 = load(sprintf('%s%.0f_ev_specific.mat', path2, sp));
    note_data2 = load(sprintf('%s%.0f_note_events.mat', path2, sp));
    state_data2 = load(sprintf('%s%.0f_state.mat', path2, sp));
    % dist_data2 = load(sprintf('%s%.0f_distances.mat', path2, sp));

    data2.(sp_str).aircraft = ac_data2;
    % data2.(sp_str).bez = bez_data2;
    data2.(sp_str).delay = delay_data2;
    % data2.(sp_str).dub = dub_data2;
    data2.(sp_str).ev_specific = ev_data2;
    data2.(sp_str).note_events = note_data2;
    data2.(sp_str).state = state_data2;
    % data2.(sp_str).distances = dist_data2;

    dist_data3 = load(sprintf('%s%.0f_distances.mat', path3, sp));
    data3.(sp_str).distances = dist_data3;
    ev_dat3 = load(sprintf('%s%.0f_ev_specific.mat', path3, sp));
    data3.(sp_str).ev_specific = ev_data2;

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
delay_matrix2 = zeros(length(xl), num_aircraft);

% Loop through each sp to populate the delay matrix
for i = 1:length(xl)
    
    sp_str = sprintf('sp%.0f', xl(i)); % Create the dynamic field name
    disp(sp_str);
    % disp(min(data.(sp_str).ev_specific.ev_specific_dat.TOI_Dist))
    min_dist = min(data.(sp_str).ev_specific.ev_specific_dat.TOI_Dist);
    % % min_dist = min(data.(sp_str).distances.dist_dat.ac_dist);
    min_dist2 = min(data2.(sp_str).ev_specific.ev_specific_dat.TOI_Dist);
    % min_dist2 = min(data2.(sp_str).distances.dist_dat.ac_dist);
    min_dist3 = min(data3.(sp_str).ev_specific.ev_specific_dat.TOI_Dist);
    % min_dist3 = min(data3.(sp_str).distances.dist_dat.ac_dist);
    if data.(sp_str).delay.delay_dat.delay >= 0
    delay_values = data.(sp_str).delay.delay_dat.delay; % Extract delay values for this sp
    delay_values2 = data2.(sp_str).delay.delay_dat.delay;
    % if xl(i) == 311
    %     delay_values = delay_values(1:2);
        % disp(delay_values)
    % end
    
    delay_matrix(i, :) = delay_values(); 
    mins(i) = min_dist;
    delay_matrix2(i, :) = delay_values2; 
    mins2(i) = min_dist2;
    mins3(i) = min_dist3;
    end

   % Store the delay values in the matrix
end
figure(1);
% Create bar plot with grouped bars
x = 1:length(xl); % X values corresponding to sp indices
% bar(x, [delay_matrix, delay_matrix2], bar_width); % Plot grouped bars
bar(x, delay_matrix, bar_width); % Plot grouped bars

% Customize the plot
xticks(x);
% grid();
xticklabels(arrayfun(@num2str, xl, 'UniformOutput', false)); % Set x-axis labels to sp numbers
xlabel('Initial Spacing Between Fleet Aircraft (m)');
ylabel('Delay (s)');
legend({'AX0', 'AX1'}, 'Location', 'best');
% title('Delay for Each Aircraft Across SP Values');
ylim([0, 25]);
grid('on');
hold off;


figure(2);
% bar(x, mins, bar_width); % Plot grouped bars% home = [39.49029047106718, -82.2];
b=bar(x, [mins; mins2;]', bar_width);
% b(1).FaceColor = 'green';
% b(2).FaceColor = 'black';
% plot(x, mins, LineWidth=5)
xticks(x);
xticklabels(arrayfun(@num2str, xl, 'UniformOutput', false)); % Set x-axis labels to sp numbers
xlabel('Initial Spacing Between Aircraft (m)');
ylabel('Minimum Distance Between Aircraft (m)');
% legend({'Holding Pattern','Opposite Sides', 'Same Side'}, 'Location', 'best', 'FontSize', 14);
% title('Distance Between Fleet Aircraft and EX0 Across SP Values');
legend({'Holding Pattern','Bezier'}, 'Location', 'best', 'FontSize', 14);
ylim([0, 25]);
% xlim([40, 8])
grid('on')
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


