clear;

% System definition

A = [1 1; 0 1];
B = [0; 1];

n = 2;
m = 1;


% Standard MPC ---------------------------------------------------------------------------

% Hyperparameters
N = 6; % horizon length N
s = 1; % number of control steps actually applied s
Q = eye(n); % state cost matrix
R = eye(m); % input cost matrix
x_0 = [0.4; 0.4]; % initial state x_0
iterations = 35; % number of algorithm iterations


% cost function: 1/2 * x'Hx + f'x, subject to: Aeq * x = beq

% Build H matrix and f vector
Q_d= zeros(N*n, N*n);
R_d = zeros(N*m, N*m);
for i = 1:N
   Q_d((i-1)*n+1:(i-1)*n+n, (i-1)*n+1:(i-1)*n+n) = Q;
   R_d((i-1)*m+1:(i-1)*m+m, (i-1)*m+1:(i-1)*m+m) = R;
end
H = 2*blkdiag(R_d, Q_d, Q);

f = zeros(N*m+(N+1)*n, 1);

% Build Aeq and beq
Aeq = zeros((N+1)*n, N*m+(N+1)*n);
for i = 1:N
    Aeq((i-1)*n+1:(i-1)*n+n, (i-1)*m+1:(i-1)*m+m) = B;
    Aeq((i-1)*n+1:(i-1)*n+n, (i-1)*n+1+N*m:(i-1)*n+1+N*m+(n-1)) = A;
    Aeq((i-1)*n+1:(i-1)*n+n, (i-1)*n+1+N*m+n:(i-1)*n+1+N*m+n+(n-1)) = -eye(n);
end
Aeq(N*n+1:(N+1)*n, N*m+1:N*m+n) = eye(n);

beq = zeros((N+1)*n, 1);

% Build constraints
ub = [ones(N*m, 1); ones((N+1)*n, 1)];
lb = [-ones(N*m, 1); -ones((N+1)*n, 1)];
%ub = [];
%lb = [];

u = [];
x = [];
x_start = x_0;

% Controller (MPC)
for i = 1:iterations
    % update beq with the starting state
    beq(N*n+1:(N+1)*n) = x_start;

    [sol, fval, exitflag, output] = quadprog(H, f, [], [], Aeq, beq, lb, ub);
    disp(i);

    u_pred = reshape(sol(1:N*m), [m,N]);
    x_pred = reshape(sol(N*m+1:end), [n,N+1]);

    u = [u u_pred(:, 1:s)];
    x = simulate_sys(A, B, u, x_0, n);

    x_start = x(:, end);

end


figure();
scatter(0:1:size(x,2)-1, x(1,:));

figure();
scatter(0:1:size(x,2)-1, x(2,:));

figure();
plot(x(1,:), x(2,:), '-o');


% Tube MPC -------------------------------------------------------------------------------
% Computation of tightened constraints

K = [-0.4 -1.2];
A_k = A+B*K;

% original constraints as scalar inequalities
c_1 = [1; 0; 0];
d_1 = 1;
c_2 = [-1; 0; 0];
d_2 = 1;
c_3 = [0; 1; 0];
d_3 = 1;
c_4 = [0; -1; 0];
d_4 = 1;
c_5 = [0; 0; 1];
d_5 = 1;
c_6 = [0; 0; -1];
d_6 = 1;

% Computation of theta_N
Ns = 25;
theta_N = [];

for N = 1:Ns
    lb = -0.1*ones(N*n, 1);
    ub = 0.1*ones(N*n, 1);

    f_lin_1 = built_f_lin_vector(c_1, K, A_k, n, N);
    f_lin_2 = built_f_lin_vector(c_2, K, A_k, n, N);
    f_lin_3 = built_f_lin_vector(c_3, K, A_k, n, N);
    f_lin_4 = built_f_lin_vector(c_4, K, A_k, n, N);
    f_lin_5 = built_f_lin_vector(c_5, K, A_k, n, N);
    f_lin_6 = built_f_lin_vector(c_6, K, A_k, n, N);

    [x1, fval1] = linprog(f_lin_1, [], [], [], [], lb, ub);
    [x2, fval2] = linprog(f_lin_2, [], [], [], [], lb, ub);
    [x3, fval3] = linprog(f_lin_3, [], [], [], [], lb, ub);
    [x4, fval4] = linprog(f_lin_4, [], [], [], [], lb, ub);
    [x5, fval5] = linprog(f_lin_5, [], [], [], [], lb, ub);
    [x6, fval6] = linprog(f_lin_6, [], [], [], [], lb, ub);

    theta_N = [theta_N [-fval1; -fval2; -fval3; -fval4; -fval5; -fval6]];
end


% Computation of alpha
w_set = [0.1 0.1 -0.1 -0.1; 0.1 -0.1 0.1 -0.1];
alpha = [];

for N = 0:Ns
    alpha_i = max([max_inf_norm(A_k^N, w_set)/max_inf_norm(eye(n), w_set) ... 
                    max_inf_norm(K*A_k^N, w_set)/max_inf_norm(K, w_set)]);

    alpha = [alpha alpha_i];
end

figure();
semilogy(0:1:Ns, alpha, '-o');
xlabel('$N$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$\alpha$', 'Interpreter', 'latex', 'FontSize', 16);

% Computation of new (tightened) bounds
% d - (1-alpha)^-1 theta_N
d = repmat([d_1; d_2; d_3; d_4; d_5; d_6], 1, Ns);
d_tightened = d - ones(6, Ns)./(ones(6, Ns) - repmat(alpha(2:end), 6, 1)).*theta_N;

figure();
plot(3:1:Ns, d_tightened(1, 3:end), '-o');
hold on
plot(3:1:Ns, d_tightened(3, 3:end), '-s');
hold on
plot(3:1:Ns, d_tightened(5, 3:end), '-^');
xlabel('$N$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('component-wise tightened bound', 'Interpreter', 'latex', 'FontSize', 16);

% So the following tightened constraints can be selected (from the preceding plot at N=16):
% |x1| =< 0.445203
% |x2| =< 0.698194
% |u| =< 0.747355


% Hyperparameters
N_t = 6; % horizon length N
s_t = 1; % number of control steps actually applied s
Q_t = eye(n); % state cost matrix
R_t = eye(m); % input cost matrix
x_0_t= [0.4; -0.6]; % initial state x_0
iterations_t = 35; % number of algorithm iterations


% cost function: 1/2 * x'Hx + f'x, subject to: Aeq * x = beq

% Build H matrix and f vector
Q_d_t= zeros(N_t*n, N_t*n);
R_d_t = zeros(N_t*m, N_t*m);
for i = 1:N_t
   Q_d_t((i-1)*n+1:(i-1)*n+n, (i-1)*n+1:(i-1)*n+n) = Q_t;
   R_d_t((i-1)*m+1:(i-1)*m+m, (i-1)*m+1:(i-1)*m+m) = R_t;
end
H_t = 2*blkdiag(R_d_t, Q_d_t, Q_t);

f_t = zeros(N_t*m+(N_t+1)*n, 1);

% Build Aeq and beq
Aeq_t = zeros((N_t+1)*n, N_t*m+(N_t+1)*n);
for i = 1:N_t
    Aeq_t((i-1)*n+1:(i-1)*n+n, (i-1)*m+1:(i-1)*m+m) = B;
    Aeq_t((i-1)*n+1:(i-1)*n+n, (i-1)*n+1+N_t*m:(i-1)*n+1+N_t*m+(n-1)) = A;
    Aeq_t((i-1)*n+1:(i-1)*n+n, (i-1)*n+1+N_t*m+n:(i-1)*n+1+N_t*m+n+(n-1)) = -eye(n);
end
Aeq_t(N_t*n+1:(N_t+1)*n, N_t*m+1:N_t*m+n) = eye(n);

beq_t = zeros((N_t+1)*n, 1);

% Build constraints
lb_t = [-ones(N_t*m, 1)*0.747355; repmat([-0.445203; -0.698194], N_t+1, 1)];
ub_t = [ones(N_t*m, 1)*0.747355; repmat([0.445203; 0.698194], N_t+1, 1)];
%ub_t = [];
%lb_t = [];

u_t = [];
x_t = [];
x_start_t = x_0_t;

% Controller (MPC) (nominal system)
for i = 1:iterations_t
    % update beq with the starting state
    beq_t(N_t*n+1:(N_t+1)*n) = x_start_t;

    [sol_t, fval_t, exitflag_t, output_t] = quadprog(H_t, f_t, [], [], Aeq_t, beq_t, lb_t, ub_t);
    disp(i);

    u_pred_t = reshape(sol_t(1:N_t*m), [m,N_t]);
    x_pred_t = reshape(sol(N_t*m+1:end), [n,N_t+1]);

    u_t = [u_t u_pred_t(:, 1:s_t)];
    x_t = simulate_sys(A, B, u_t, x_0_t, n);

    x_start_t = x_t(:, end);

end


figure();
scatter(0:1:size(x_t,2)-1, x_t(1,:));
xlabel('time step', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$\bar{x}_{1}$', 'Interpreter', 'latex', 'FontSize', 16);
grid on;

figure();
scatter(0:1:size(x_t,2)-1, x_t(2,:));
xlabel('time step', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$\bar{x}_{2}$', 'Interpreter', 'latex', 'FontSize', 16);
grid on;

figure();
plot([0 0], [1 -1], '-', 'Color', [0.5 0.5 0.5]);
hold on
plot([1 -1], [0 0], '-', 'Color', [0.5 0.5 0.5]);
hold on
plot(x_t(1,:), x_t(2,:), '-o', 'Color', [0 0.4470 0.7410]);
hold on
plot(x_t(1,1), x_t(2,1), '.', 'Color', [0 0.4470 0.7410], 'MarkerSize', 30);
ax = gca;
ax.FontSize = 16;
xlabel('\boldmath$\bar{x}_{1}$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('\boldmath$\bar{x}_{2}$', 'Interpreter', 'latex', 'FontSize', 20, 'Rotation', 0);
xlim([-0.5 0.5]);
ylim([-0.7 0.3]);
%grid on;

% Simulate uncertain system (with additive bounded disturbance)
w_bound = 0.1;
[x_uncertain, u_uncertain] = simulate_uncertain_sys(A, B, u_t, x_t, x_0_t, n, m, K, w_bound);

figure();
scatter(0:1:size(x_uncertain,2)-1, x_uncertain(1,:));
xlabel('time step', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$x_{1}$', 'Interpreter', 'latex', 'FontSize', 16);
grid on;

figure();
scatter(0:1:size(x_uncertain,2)-1, x_uncertain(2,:));
xlabel('time step', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$x_{2}$', 'Interpreter', 'latex', 'FontSize', 16);
grid on;

figure();
plot([0 0], [1 -1], '-', 'Color', [0.5 0.5 0.5]);
hold on
plot([1 -1], [0 0], '-', 'Color', [0.5 0.5 0.5]);
hold on
plot(x_uncertain(1,:), x_uncertain(2,:), '-o', 'Color', [0 0.4470 0.7410]);
hold on
plot(x_uncertain(1,1), x_uncertain(2,1), '.', 'Color', [0 0.4470 0.7410], 'MarkerSize', 30);
ax = gca;
ax.FontSize = 16;
xlabel('\boldmath$x_{1}$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('\boldmath$x_{2}$', 'Interpreter', 'latex', 'FontSize', 20, 'Rotation', 0);
xlim([-0.5 0.5]);
ylim([-0.7 0.3]);
%grid on;

% Function definitions
function x = simulate_sys(A, B, u, x_0, n)
    x = [x_0 zeros(n, size(u,2))];
    for i = 2:1+size(u,2)
        x(:, i) = A*x(:,i-1) + B*u(:,i-1);
    end
end

function f_lin = built_f_lin_vector(c, K, A_k, n, N)
    f_lin = zeros(N*n, 1);
    for i = 0:N-1
        f_lin(i*n+1:i*n+n, 1) = -(c.' * [eye(n); K] * A_k^(i)).';
    end
end

function sol = max_inf_norm(F, w_set)
    sol = zeros(1, size(w_set,2));
    for i=1:size(w_set, 2)
        sol(1, i) = norm(F*w_set(:,i), Inf);
    end
    sol = max(sol);
end

function [x, u] = simulate_uncertain_sys(A, B, u_n, x_n, x_0, n, m, K, w_bound)
    x = [x_0 zeros(n, size(u_n,2))];
    u = zeros(m, size(u_n,2));
    for i = 2:1+size(u_n,2)
        u(:,i-1) = u_n(:,i-1) + K*(x(:,i-1)-x_n(:,i-1));
        x(:,i) = A*x(:,i-1) + B*u(:,i-1) + unifrnd(-w_bound, w_bound, n, 1);
    end
end

