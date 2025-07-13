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

% Initialize the block Laplacian matrix (size: 20x20)
L = zeros(n * d);

% Helper function to compute (e_i - e_j)(e_i - e_j)^T ⊗ W
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
        0.5 * cos(0.3 * t),       0.4 * cos(0.2 * i * t + 0.15)];
end

Q = @(t) buildQ(Qblocks, t);

%% first_order
T = 0.1; alpha2 = 10 ;alpha1 = 14.95 ;                          % sampling_time, gains of control input
x = [36;13;14;2;9;-42;26;-37;-31;-26;-8;-45;41;45;0;-1;-16;41;-13;-38] ;                    % x(1) : state
v = randi(10,20,1) -5 ; 
z = zeros(20,1) ;                            %z(1) : estimation of state
d_hat= zeros(20,1)  ;
d= [0.15;1;0.05;0.8;0.85;0.9;0.1;0.4;0.3;0.85;0.45;0.95;0.2;0.3;0.15;0.15;0.9;0.6;0.55;0.15];

%record the state
x_record = x; v_record = v;
d_hat_record = d_hat;
z = d_hat - pinv(Q(T))*(-alpha1 *T*L*x + (1- alpha2*T)*v) ;
for i = 1: 1000


    [x_new, v_new] = update_state_second_order(x,v, alpha1, alpha2, T, L, Q(i*T), d, d_hat );
    x= x_new; x_record = [x_record; x];
    v= v_new; v_record = [v_record; v];

    [d_hat, z_new] = disturbance_observer_second_order( z, pinv(Q(T*i)),x, v, pinv(Q((i+1)*T)), L, alpha1, alpha2, T) ;
    d_hat_record =[d_hat_record;d_hat];
    z=z_new;
    

end

% plot results
figure(1); clf;
hold on;

figure(1); clf;
hold on;

for i = 1:2:19
    plot(x_record(i:20:20000), 'k'); 
end
for i = 2:2:20
    plot(x_record(i:20:20000), 'b');  
end

h1 = plot(NaN, NaN, 'k', 'LineWidth', 1.5);
h2 = plot(NaN, NaN, 'b', 'LineWidth', 1.5);

lgd = legend([h1 h2], ...
    {'$x_{i,1},\ i=1,\ldots,10$', '$x_{i,2},\ i=1,\ldots,10$'}, ...
    'Interpreter','latex', ...
    'Location','northeast');

lgd.FontSize = 20;
lgd.ItemTokenSize = [30, 10];  

xlabel('$k\ [\mathrm{iteration}]$', 'Interpreter','latex','FontSize',18);
ylabel('$x_{i,k},\ i=1,\ldots,10,\ k=1,2$', 'Interpreter','latex','FontSize',18);
axis tight;
grid on;

exportgraphics(gcf, 'second_order_matched_state3.pdf', ...
    'BackgroundColor', 'white', ...
    'ContentType', 'vector'); 


figure(2); clf;
hold on;

for i = 1:2:19
    plot(v_record(i:20:10000), 'k');  
end
for i = 2:2:20
    plot(v_record(i:20:10000), 'b'); 
end

h1 = plot(NaN, NaN, 'k', 'LineWidth', 1.5);
h2 = plot(NaN, NaN, 'b', 'LineWidth', 1.5);

lgd = legend([h1 h2], ...
    {'$v_{i,1},\ i=1,\ldots,10$', '$v_{i,2},\ i=1,\ldots,10$'}, ...
    'Interpreter','latex', ...
    'Location','northeast');

lgd.FontSize = 20;
lgd.ItemTokenSize = [30, 10];

xlabel('$k\ [\mathrm{iteration}]$', 'Interpreter','latex','FontSize',18);
ylabel('$v_{i,k},\ i=1,\ldots,10,\ k=1,2$', 'Interpreter','latex','FontSize',18);
axis tight;
grid on;

exportgraphics(gcf, 'second_order_matched_velocity3.pdf', ...
    'BackgroundColor', 'white', ...
    'ContentType', 'vector'); 



figure(3); clf;
hold on;

for i = 1:2:19
    plot(d_hat_record(i:20:i+20*3), 'b');  
end

for i = 1:2:length(d)
    yline(d(i), 'k--');  
end

h1 = plot(NaN, NaN, 'b-');   
h2 = plot(NaN, NaN, 'k--');    

lgd = legend([h1 h2], ...
    {'$\hat{d}_{i,1},\ i=1,\ldots,10$', '$d_{i,1},\ i=1,\ldots,10$'}, ...
    'Interpreter','latex', ...
    'Location','northeast');

lgd.FontSize = 20;
lgd.ItemTokenSize = [30, 10];

xlabel('$k\ [\mathrm{iteration}]$', 'Interpreter','latex','FontSize',18);
ylabel('$\hat{d}_{i,1},\ d_{i,1},\ i=1,\ldots,10$', 'Interpreter','latex','FontSize',18);

ylim([0 1.3]);  
xticks([1 2 3]);
axis normal;
grid on;
exportgraphics(gcf, 'second_order_matched_first_disturbance3.pdf', ...
    'BackgroundColor', 'white', ...
    'ContentType', 'vector');

figure(4); clf;
hold on;

for i = 2:2:20
    plot(d_hat_record(i:20:i+20*3), 'b'); 
end

for i = 2:2:length(d)
    yline(d(i), 'k--');  end

h1 = plot(NaN, NaN, 'b-');    
h2 = plot(NaN, NaN, 'k--');   

lgd = legend([h1 h2], ...
    {'$\hat{d}_{i,2},\ i=1,\ldots,10$', '$d_{i,2},\ i=1,\ldots,10$'}, ...
    'Interpreter','latex', ...
    'Location','northeast');

lgd.FontSize = 20;
lgd.ItemTokenSize = [30, 10];

% Ghi nhãn trục
xlabel('$k\ [\mathrm{iteration}]$', 'Interpreter','latex','FontSize',18);
ylabel('$\hat{d}_{i,2},\ d_{i,2},\ i=1,\ldots,10$', 'Interpreter','latex','FontSize',18);
xticks([1 2 3]);
ylim([0 1.5]);      
axis normal;
grid on;
exportgraphics(gcf, 'second_order_matched_second_disturbance3.pdf', ...
    'BackgroundColor', 'white', ...
    'ContentType', 'vector');


figure(5);
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

lgd=legend([h_list; h_start; h_end], ...
    [labels, {'initial position', 'final position'}], ...
    'Interpreter','latex', ...
    'Location','eastoutside');

xlabel('$x_{i,1}$', 'Interpreter','latex','FontSize',18);
ylabel('$x_{i,2}$', 'Interpreter','latex','FontSize',18);

exportgraphics(gcf, 'second_order_matched_trajectories3.pdf', ...
    'BackgroundColor', 'white', ...
    'ContentType', 'vector');

%% function for system
function [d_hat , z_new] = disturbance_observer_second_order( z, L_mat_obs_old, x, v, L_mat_obs, Laplacian_mat, alpha1, alpha2, T)
d_hat = z + L_mat_obs_old * v; 
z_new = d_hat - L_mat_obs * ((1 - alpha2 * T ) *v - alpha1*T*Laplacian_mat*x);
end 

function [x_new, v_new] = update_state_second_order(x_old,v_old , alpha1, alpha2, T, Laplacian_mat, Q_matrix, d_old, d_hat_old )
    x_new = x_old + T * v_old   ;
    v_new = v_old + T * (-alpha1 * Laplacian_mat *x_old - alpha2*v_old) + Q_matrix*(d_old - d_hat_old);
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

