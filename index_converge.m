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
index_convergence = zeros(19,1);
%% first_order
eigL= eig(L);
for (vonglap = 1:19)
T = 0.1; alpha = (0.05*vonglap) * 2/(T*eigL(end));                           % sampling_time, gain of control input
x = [36;13;14;2;9;-42;26;-37;-31;-26;-8;-45;41;45;0;-1;-16;41;-13;-38] ;                 % x(1) : state
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
for i = 1: 5000 
    
    x = update_state_first_order(x,alpha,T, L, Q(i*T), d, d_hat,u);  %calculate x(k+1)
    x_record = [x_record; x];
   
    [d_hat, z] = disturbance_observer_first_order(d_hat, z, pinv(Q((i)*T)), x, pinv(Q((i+1)*T)), L, alpha, T); %calculate d(k+1) and z(k+2)
    d_hat_record = [d_hat_record; d_hat];
     
end

index_convergence(vonglap) = find_index(x_record);
end
figure(1); clf;

alpha_vals = 0.05:0.05:0.95;

plot(alpha_vals, index_convergence, '-o', ...
    'LineWidth', 2, ...
    'MarkerSize', 6);

xlabel('$\alpha$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('Number of iterations', 'FontSize', 20);

grid on;
set(gca, 'FontSize', 12);  % chỉnh số trên trục

xlim([0 1]);
ylim([0, max(index_convergence)*1.05]);

% Xuất file vector chất lượng cao
exportgraphics(gcf, 'convergence_time.pdf', ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');


%% function for system
function [d , z_new] = disturbance_observer_first_order(d, z, L_mat_obs, x, L_mat_obs_new, L, alpha, T)
 d = z + L_mat_obs * x;
 z_new = z + L_mat_obs * x - L_mat_obs_new * (eye(size(L,1)) - alpha * T * L) *x ;
end 

function x = update_state_first_order(x,alpha,T, Laplacian_mat, Q_matrix, d, d_hat,u )
    x = (eye(size(Laplacian_mat,1)) - alpha*T*Laplacian_mat)* x + Q_matrix*(d-d_hat) - u;
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