
clear; clc;

% ==== 1. Định nghĩa các cạnh theo nhóm màu ====
black_edges = [1 2; 3 4; 4 5; 6 7; 7 8; 9 10];
red_edges   = [2 3; 5 6; 8 9];
blue_edges  = [1 3; 2 4; 4 6; 5 7; 7 9; 8 10];

all_edges = [black_edges; red_edges; blue_edges];
G = graph(all_edges(:,1), all_edges(:,2));

% ==== 2. Định nghĩa vị trí node ====
x = [-3, -2, -2, -1,  0,  0, 1,  2,  2, 3]; 
y = [ 0,  1, -1,  0,  1, -1, 0,  1, -1, 0];
nodePos = [x' y'];

% ==== 3. Vẽ đồ thị ====
figure('Color','white', 'Units','inches', 'Position',[1 1 6 2.5]); % chiều ngang dàn đều hơn
h = plot(G, ...
    'XData', nodePos(:,1), ...
    'YData', nodePos(:,2), ...
    'LineWidth', 3, ...
    'MarkerSize', 20, ...
    'NodeColor', [0.3 0.8 0.3], ...    % xanh lá cây
    'EdgeColor', [0 0 0], ...    %  màu đen 
    'NodeLabel', []);


for i = 1:numnodes(G)
    text(x(i), y(i), num2str(i), ...
        'FontSize', 14, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle');
end

highlight(h, red_edges(:,1), red_edges(:,2), 'EdgeColor', [1 0 0],'LineWidth',3);
highlight(h, blue_edges(:,1), blue_edges(:,2), 'EdgeColor', [0.2 0.5 1.0],'LineWidth',3);


axis equal off;
set(gca, 'Position', [0 0 1 ]);
% ==== 7. Lưu hình ảnh ====
exportgraphics(gcf, 'graph.pdf', ...
    'BackgroundColor', 'white', ...
    'ContentType', 'vector'); 


