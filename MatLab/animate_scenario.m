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

AC0_lon2 = AC0(1, 1:1:401);
AC0_lat2 = AC0(2, 1:1:401);
AC1_lon2 = AC1(1, 1:1:401);
AC1_lat2 = AC1(2, 1:1:401);
EM0_lon2 = EM0(1, 1:1:401);
EM0_lat2 = EM0(2, 1:1:401);


% subplot(3, 3, 1);
rect_handle = []; 
AC0_text = [];
AC1_text = [];
EM0_text = [];
AC0_q = [];
AC1_q = [];
EM0_1 = [];

box('on')
ax = gca;
ax.FontWeight = 'bold';
ax.FontName = 'Times';
set(gcf, 'Position', [700 100 1000 1100])
axis('equal');
grid();
xlim([EM0_lon(1)-0.005, EM0_lon(401)+0.005])
% ylim([AC1_lat(1)-0.0015, EM0_lat(401)+0.0001])
ylim([39.417, 39.4425])
ylabel('Latitude', FontName='Times')
xlabel('Longitude', FontName='Times')
r = 0.0004144027532220207;
num_frames = size(AC0_lat2, 2);
frames(2:num_frames) = struct('cdata', [], 'colormap', []);
for i = 2:1:num_frames
    % pause(0.005)
    delete(rect_handle);
    delete(AC0_text);
    delete(AC1_text);
    delete(EM0_text);
    delete(AC0_q);
    delete(AC1_q);
    delete(EM0_1);

    hold on;
    plot([AC0_lon2(i-1), AC0_lon2(i)], [AC0_lat2(i-1), AC0_lat2(i)], LineWidth=4, color = 'blue')
    plot([AC1_lon2(i-1), AC1_lon2(i)], [AC1_lat2(i-1), AC1_lat2(i)], LineWidth=4, color = 'green')
    plot([EM0_lon2(i-1), EM0_lon2(i)], [EM0_lat2(i-1), EM0_lat2(i)], LineWidth=4, color = 'red')

    scatter_handles = findobj(gca, 'Type', 'Scatter'); % Find existing scatter plots
    delete(scatter_handles);


    AC0_q = quiver(AC0_lon(i), AC0_lat(i), cosd(90-AC0_hdg(i)), sind(90-AC0_hdg(i)), 0.001, 'Color', 'black', 'LineWidth', 2, 'MaxHeadSize',2);
    AC1_q = quiver(AC1_lon(i), AC1_lat(i), cosd(90-AC1_hdg(i)), sind(90-AC1_hdg(i)), 0.001, 'Color', 'black', 'LineWidth', 2, 'MaxHeadSize',2);
    EM0_1 = quiver(EM0_lon(i), EM0_lat(i), 0, 1, 0.001, 'Color', 'black', 'LineWidth', 2, 'MaxHeadSize',2);

    scatter(AC0_lon(i), AC0_lat(i), 50, 'o', 'filled', 'blue')
    scatter(EM0_lon(i), EM0_lat(i), 50, '^',  'filled', 'red')
    scatter(AC1_lon(i), AC1_lat(i), 50, 'o', 'filled', 'green')
    rect_handle = rectangle('Position', [EM0_lon(i)-r, EM0_lat(i)-r, 2*r, 2*r], 'Curvature', [1, 1], 'EdgeColor','r', 'LineWidth',2);

    AC0_text = text(AC0_lon(i)+1*r, AC0_lat(i), sprintf('AX0'), 'FontSize', 18, FontName='Times', FontWeight= 'bold');
    AC1_text = text(AC1_lon(i)+1*r, AC1_lat(i), sprintf('AX1'), 'FontSize', 18, FontName='Times', FontWeight= 'bold');
    EM0_text = text(EM0_lon(i)-5*r, EM0_lat(i), sprintf('EX0'), 'FontSize', 18, FontName='Times', FontWeight= 'bold');

    
    drawnow;
    frames(i) = getframe(gcf);
    hold off;
    % hold on;
    
    

end
% % Convert frames to indexed images
% [im, map] = rgb2ind(frames(1).cdata, 256);
% im(1, 1, 1, num_frames) = 0;
% for i = 1:num_frames
%     [im(:,:,1,i), map] = rgb2ind(frames(i).cdata, map);
% end
% 
% % Save frames as GIF
% filename = 'Scen1.gif';
% for i = 1:num_frames
%     if i == 1
%         imwrite(im(:,:,1,i), map, filename, 'gif', 'Loopcount', inf, 'DelayTime', .5);
%     else
%         imwrite(im(:,:,1,i), map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', .0000001);
%     end
% end
% Convert frames to indexed images
frames_to_convert = frames(2:num_frames); % Only use valid frames
[im, map] = rgb2ind(frames_to_convert(1).cdata, 256);
im(1, 1, 1, numel(frames_to_convert)) = 0; % Preallocate GIF array
for k = 1:numel(frames_to_convert)
    im(:, :, 1, k) = rgb2ind(frames_to_convert(k).cdata, map);
end

% Save frames as GIF
filename = 'Scen1.gif';
for k = 1:numel(frames_to_convert)
    if k == 1
        imwrite(im(:, :, 1, k), map, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.005);
    else
        imwrite(im(:, :, 1, k), map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.005);
    end
end