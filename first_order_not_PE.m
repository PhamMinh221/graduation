clear;clc;
code = 0;

clear;clc;
    black_edges = [1 2; 3 4; 4 5; 6 7; 7 8; 9 10];
    red_edges = [2 3; 5 6; 8 9];
    blue_edges = [1 3; 2 4; 4 6; 5 7; 7 9; 8 10];
% Number of nodes and block size
n = 10;      % total number of nodes
d = 2;       % dimension of each node (block size)

% Matrix weights for different edge types (each 2x2)
W_black = 0.5 * [1 -1; -1 1];
W_red   = 0.5 * eye(2);
W_blue  = 0.5 * [1 1; 1 1];


L = zeros(n * d);

add_edge_block = @(i, j, W) ...
    kron((unit_vector(i,n) - unit_vector(j,n)) * ...
         (unit_vector(i,n) - unit_vector(j,n))', W);

% Add contributions from black edges
for k = 1:size(black_edges,1)
    i = black_edges(k,1);
    j = black_edges(k,2);
    L = L + add_edge_block(i, j, W_black);
end

% Add contributions from red edges
for k = 1:size(red_edges,1)
    i = red_edges(k,1);
    j = red_edges(k,2);
    L = L + add_edge_block(i, j, W_red);
end

% Add contributions from blue edges
for k = 1:size(blue_edges,1)
    i = blue_edges(k,1);
    j = blue_edges(k,2);
    L = L + add_edge_block(i, j, W_blue);
end

Qblocks = cell(n, 1);
for i = 1:n
    Qblocks{i} = @(t) [ ...
        sin(0.1 * t * i),         0.2 * exp(-i * t);
        3 * sin(0.1 * t * i),       0.6 * exp(-i * t)];
end

Q = @(t) buildQ(Qblocks, t);

%% first_order
eigL= eig(L);

T = 0.1; alpha = 0.3 * 2/(T*eigL(end));                           % sampling_time, gain of control input
x = [36;13;14;2;9;-42;26;-37;-31;-26;-8;-45;41;45;0;-1;-16;41;-13;-38] ;                  % x(1) : state
z = zeros(20,1) ;                            %z(1) : estimation of state
d_hat= ones(20,1)  ;
d=  [0.15;1;0.05;0.8;0.85;0.9;0.1;0.4;0.3;0.85;0.45;0.95;0.2;0.3;0.15;0.15;0.9;0.6;0.55;0.15];
% update 
[d_hat, z] = disturbance_observer_first_order(d_hat, z, zeros(size(L,1)), x,pinv(Q(T)), L, alpha, T); % calcultate d(1) and z(2)
z_check = z;
%record the state
x_record = x;
d_hat_record = d_hat;
u=0;
u_new =0;
for i = 1: 500    
    x = update_state_first_order(x,alpha,T, L, Q(i*T), d, d_hat);  %calculate x(k+1)
    x_record = [x_record; x];
     if i < 3
    [d_hat, z] = disturbance_observer_first_order(d_hat, z, pinv(Q((i)*T)), x, pinv(Q((i+1)*T)), L, alpha, T); %calculate d(k+1) and z(k+2)
    d_hat_record = [d_hat_record; d_hat];
     end 
end

index_convergence = find_index(x_record);

% plot results
figure(1); clf;
hold on;


for i = 1:2:19
    plot(x_record(i:20:10000), 'k');  
end
for i = 2:2:20
    plot(x_record(i:20:10000), 'b'); 
end

h1 = plot(NaN, NaN, 'k', 'LineWidth', 1.5);
h2 = plot(NaN, NaN, 'b', 'LineWidth', 1.5);


lgd1 = legend([h1 h2], ...
    {'$x_{i,1},\ i=1,\ldots,10$', '$x_{i,2},\ i=1,\ldots,10$'}, ...
    'Interpreter','latex', ...
    'Location','northeast');
lgd1.FontSize = 20;
lgd1.ItemTokenSize = [30, 10];

xlabel('$k\ [\mathrm{iteration}]$', 'Interpreter','latex', 'FontSize', 18);
ylabel('$x_{i,k},\ i=1,\ldots,10,\ k=1,2$', 'Interpreter','latex', 'FontSize', 18);

set(gca, 'FontSize', 16);

axis tight;
grid on;

exportgraphics(gcf, 'first_order_state3.pdf', ...
    'BackgroundColor', 'white', ...
    'ContentType', 'vector');


figure(2); clf;
hold on;
for i = 1:2:19
    plot(1:3, d_hat_record(i:20:end), 'b');  
end

for i = 1:2:length(d)
    yline(d(i), 'k--');  
end

h1 = plot(NaN, NaN, 'b-');     
h2 = plot(NaN, NaN, 'k--');  

lgd2 = legend([h1 h2], ...
    {'$\hat{d}_{i,1},\ i=1,\ldots,10$', '$d_{i,1},\ i=1,\ldots,10$'}, ...
    'Interpreter','latex', ...
    'Location','northeast');

lgd2.FontSize = 20;
lgd2.ItemTokenSize = [30, 10];

xlabel('$k\ [\mathrm{iteration}]$', 'Interpreter','latex','FontSize',18);
ylabel('$\hat{d}_{i,1},\ d_{i,1},\ i=1,\ldots,10$', 'Interpreter','latex','FontSize',18);
set(gca, 'FontSize', 16); 
ylim([0 1.3]);      
xticks([1 2 3]);
axis normal;
grid on;
exportgraphics(gcf, 'first_order_first_disturbance3.pdf', ...
    'BackgroundColor', 'white', ...
    'ContentType', 'vector');

figure(3); clf;
hold on;

for i = 2:2:20
    plot(d_hat_record(i:20:end), 'b');  
end

for i = 2:2:length(d)
    yline(d(i), 'k--');  
end

h1 = plot(NaN, NaN, 'b-');     
h2 = plot(NaN, NaN, 'k--');    

lgd3 = legend([h1 h2], ...
    {'$\hat{d}_{i,2},\ i=1,\ldots,10$', '$d_{i,2},\ i=1,\ldots,10$'}, ...
    'Interpreter','latex', ...
    'Location','northeast');

lgd3.FontSize = 20;
lgd3.ItemTokenSize = [30, 10];

% Ghi nhãn trục
xlabel('$k\ [\mathrm{iteration}]$', 'Interpreter','latex','FontSize',18);
ylabel('$\hat{d}_{i,2},\ d_{i,2},\ i=1,\ldots,10$', 'Interpreter','latex','FontSize',18);
set(gca, 'FontSize', 16); 
ylim([0 1.5]);
xticks([1 2 3]);
axis normal;
grid on;
exportgraphics(gcf, 'first_order_second_disturbance3.pdf', ...
    'BackgroundColor', 'white', ...
    'ContentType', 'vector');


figure(4); clf;
hold on;
grid on;
axis equal;

colors = lines(10);     
line_width = 2.5;

h_list = gobjects(10,1);  

for i = 1:10
    x = x_record(2*i-1:20:end);  
    y = x_record(2*i:20:end);    

    h_list(i) = plot(x, y, '-', ...
        'Color', colors(i,:), 'LineWidth', line_width);

    plot(x(1), y(1), 'o', ...
        'MarkerSize', 10, 'LineWidth', 1.5, ...
        'Color', colors(i,:), 'MarkerFaceColor', 'none');

    plot(x(end), y(end), 's', ...
        'MarkerSize', 10, 'LineWidth', 1.5, ...
        'Color', colors(i,:), 'MarkerFaceColor', 'none');
end

h_start = plot(NaN, NaN, 'ko', ...
    'MarkerSize', 10, 'LineWidth', 1.5, 'MarkerFaceColor', 'none');
h_end   = plot(NaN, NaN, 'ks', ...
    'MarkerSize', 10, 'LineWidth', 1.5, 'MarkerFaceColor', 'none');

labels = arrayfun(@(i) sprintf('$x_{%d}$', i), 1:10, 'UniformOutput', false);

legend([h_list; h_start; h_end], ...
    [labels, {'initial position', 'final position'}], ...
    'Interpreter','latex', ...
    'Location','eastoutside');

xlabel('$x_{i,1}$', 'Interpreter','latex');
ylabel('$x_{i,2}$', 'Interpreter','latex');

exportgraphics(gcf, 'trajectories_first_order.pdf', ...
    'BackgroundColor', 'white', ...
    'ContentType', 'vector');



%% function for system
function [d , z_new] = disturbance_observer_first_order(d, z, L_mat_obs, x, L_mat_obs_new, L, alpha, T)
 d = z + L_mat_obs * x;
 z_new = z + L_mat_obs * x - L_mat_obs_new * (eye(size(L,1)) - alpha * T * L) *x ;
end 

function x = update_state_first_order(x,alpha,T, Laplacian_mat, Q_matrix, d, d_hat )
    x = (eye(size(Laplacian_mat,1)) - alpha*T*Laplacian_mat)* x + Q_matrix*(d-d_hat) ;
end

function e = unit_vector(i, n)
    % Create the i-th standard basis vector in R^n
    e = zeros(n,1);
    e(i) = 1;
end

function Qt = buildQ(Qblocks, t)
    blocks = cellfun(@(f) f(t), Qblocks, 'UniformOutput', false);  
    Qt = blkdiag(blocks{:});  
end

function vonglap = find_index(x_record)
    so_tac_tu = 10;
    so_trang_thai = 2;
    so_phantu_moi_vong = so_tac_tu * so_trang_thai;

    tong_so_vonglap = length(x_record) / so_phantu_moi_vong;

    vonglap = -1;

    for k = 1:tong_so_vonglap
        % Trích trạng thái tại vòng lặp k
        idx_start = (k - 1) * so_phantu_moi_vong + 1;
        idx_end = idx_start + so_phantu_moi_vong - 1;
        snapshot = x_record(idx_start:idx_end);
        x_matrix = reshape(snapshot, [so_trang_thai, so_tac_tu])';
        
        trungbinh = mean(x_matrix, 1);        
        max_lech = max(abs(x_matrix - trungbinh), [], 1);

        if all(max_lech < 0.01 * abs(trungbinh))
            vonglap = k;
            return;
        end
    end
end
