clc; clear; close all;
load 311_aircraft.mat;
load 311_bez.mat;
load 311_dub.mat;
load 311_state.mat;

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










home = [39.49029047106718, -82.2];

[cx1, cy1] = quadBezier([0, 574.12], [205.6159, 564.0093], [221, 711.62], 100);

[bezlon1, bezlat1] = Meters_To_WSG84(cx1, cy1, home)


lon = state_dat.lon(13550:1:17001);
lat = state_dat.lat(13550:1:17001);



hold on;
plot(lon, lat, LineWidth=2);
plot(bezlon1, bezlat1, LineWidth=2)
xlim([-82.2-0.003, -82.2+0.003])
axis('equal')


