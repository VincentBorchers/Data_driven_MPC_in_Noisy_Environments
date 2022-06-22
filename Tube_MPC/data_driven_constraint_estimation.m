clear;

% set(groot,'defaulttextinterpreter','none');  
% set(groot, 'defaultAxesTickLabelInterpreter','none');  
% set(groot, 'defaultLegendInterpreter','none');


% System definition--------------------------------------------------------

A = [1 1; 0 1];
B = [0; 1];
C = eye(2);
D = [0; 0];

n = 2;
m = 1;
p = 2;


% nominal system: x_n+ = Ax_n + Bu_n
% uncertain system : x+ = Ax + B(u_n + K(x-x_n)) + w

% Define K
K  = [-0.4 -1.2];
% Define initial state
x_0 = [0; 0];
% Define long input sequence u_n
t = 1000000;
%t = 100;
%t = 10000000;
u_bound = 1;
u_n = unifrnd(-u_bound,u_bound,m,t);
% Define an additive bounded disturbance
w_bound = 0.1;
% Define disturbance sequence
w = unifrnd(-w_bound, w_bound, n, t);

% simulate both systems
x_n = simulate_sys(A, B, u_n, x_0, n);
x = simulate_uncertain_sys(A, B, u_n, x_n, x_0, n, m, K, w);

% Compute e = x-x_n
error = x.'-x_n.';

max_1 = max(error(:,1));
max_2 = max(error(:,2));
min_1 = min(error(:,1));
min_2 = min(error(:,2));
rectangle = [[min_1 max_1 max_1 min_1 min_1].' [min_2 min_2 max_2 max_2 min_2].'];

% and plot
figure();
scatter(error(:,1), error(:,2), 'x');
grid on;
ax = gca;
ax.FontSize = 16;
xlabel('\boldmath$e_1$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('\boldmath$e_2$', 'Interpreter', 'latex', 'FontSize', 18, 'Rotation', 0);

% this gives rough approximation of S(infty)
figure();
scatter(error(:,1), error(:,2), 'x');
hold on;
plot(rectangle(:,1), rectangle(:,2), 'k--');
grid on;
ax = gca;
ax.FontSize = 16;
xlabel('\boldmath$e_1$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('\boldmath$e_2$', 'Interpreter', 'latex', 'FontSize', 18, 'Rotation', 0);
%[k_3, area_3] = convhull(rectangle,'Simplify',true);

% convex hull of error points gives better approximation of S(infty)
[k, area] = convhull(error, 'Simplify', true);
figure();
scatter(error(:,1), error(:,2), 'x');
hold on;
plot(error(k,1), error(k,2), 'k--');
grid on;
ax = gca;
ax.FontSize = 16;
xlabel('\boldmath$e_1$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('\boldmath$e_2$', 'Interpreter', 'latex', 'FontSize', 18, 'Rotation', 0);

% formula for increasing area of polygon
% xt = A*x+(1-A)*mean(x(1:end-1))
% yt = A*y+(1-A)*mean(y(1:end-1))

% increase area of rectangle by factor c^2
c = 1.1;
rectangle_increased = zeros(size(rectangle,1),2);
rectangle_increased(:,1) = c*rectangle(:,1) + (1-c)*mean(rectangle(1:end-1,1));
rectangle_increased(:,2) = c*rectangle(:,2) + (1-c)*mean(rectangle(1:end-1,2));
%[k_1, area_1] = convhull(rectangle_increased,'Simplify',true);

% increase area of convex hull by factor c^2
err_incr = zeros(size(k,1),2);
err_incr(:,1) = c*error(k,1) + (1-c)*mean(error(k(1:end-1),1));
err_incr(:,2) = c*error(k,2) + (1-c)*mean(error(k(1:end-1),2));
%[k_2, area_2] = convhull(err_incr,'Simplify',true);

% plot increased rectangle
figure();
scatter(error(:,1), error(:,2), 'x');
hold on;
plot(rectangle(:,1), rectangle(:,2), 'k--');
hold on;
plot(rectangle_increased(:,1), rectangle_increased(:,2), 'k-');
grid on;
ax = gca;
ax.FontSize = 16;
xlabel('\boldmath$e_1$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('\boldmath$e_2$', 'Interpreter', 'latex', 'FontSize', 18, 'Rotation', 0);

% plot increased convex hull
figure();
scatter(error(:,1), error(:,2), 'x');
hold on;
plot(error(k,1), error(k,2), 'k--');
hold on;
plot(err_incr(:,1), err_incr(:,2), 'k-');
grid on;
ax = gca;
ax.FontSize = 16;
xlabel('\boldmath$e_1$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('\boldmath$e_2$', 'Interpreter', 'latex', 'FontSize', 20, 'Rotation', 0);

% Make each set a Polyhedron object in MPT3
P_rect = Polyhedron(rectangle);
P_rect.irredundantVRep;
P_rect.irredundantHRep;

P_rect_incr = Polyhedron(rectangle_increased);
P_rect_incr.irredundantVRep;
P_rect_incr.irredundantHRep;

P_polygon = Polyhedron(error(k,:));
P_polygon.irredundantVRep;
P_polygon.irredundantHRep;

P_polygon_increased = Polyhedron(err_incr);
P_polygon_increased.irredundantVRep;
P_polygon_increased.irredundantHRep;
A_s = P_polygon_increased.A;
b_s = P_polygon_increased.b;

% Define original state set X
% |x|_infty <= 1
x_set = [-1 -1; 1 -1; 1 1; -1 1; -1 -1];
P_x_set = Polyhedron(x_set);

% Compute tightened state set \bar{X}
% \bar{X} = X - S(infty)
P_x_set_tight = minus(P_x_set, P_polygon_increased);
P_x_set_tight.irredundantVRep;
P_x_set_tight.irredundantHRep;

% Compute 'rough' tightened state set \bar{X}
P_x_set_tight_rough = minus(P_x_set, P_rect_incr);
P_x_set_tight_rough.irredundantVRep;
P_x_set_tight_rough.irredundantHRep;

% Plot sets
figure();
fill(x_set(:,1), x_set(:,2), [0 0.4470 0.7410], 'FaceAlpha', 0.4);
xlim([-1.3 1.3]);
ylim([-1.3 1.3]);
grid on;
axis square;
text(x_set(1,1), x_set(1,2), ['(' num2str(x_set(1,1)) ',' num2str(x_set(1,2)) ')'], 'HorizontalAlignment','right', 'FontSize', 16);
text(x_set(2,1), x_set(2,2), ['(' num2str(x_set(2,1)) ',' num2str(x_set(2,2)) ')'], 'HorizontalAlignment','left', 'FontSize', 16);
text(x_set(3,1), x_set(3,2), ['(' num2str(x_set(3,1)) ',' num2str(x_set(3,2)) ')'], 'HorizontalAlignment','left', 'FontSize', 16);
text(x_set(4,1), x_set(4,2), ['(' num2str(x_set(4,1)) ',' num2str(x_set(4,2)) ')'], 'HorizontalAlignment','right', 'FontSize', 16);
ax = gca;
ax.FontSize = 16;
xlabel('\boldmath$x_1$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('\boldmath$x_2$', 'Interpreter', 'latex', 'FontSize', 20, 'Rotation', 0);

figure();
fill(P_polygon_increased.V(:,1), P_polygon_increased.V(:,2), [0 0.4470 0.7410], 'FaceAlpha', 0.4);
xlim([-1.3 1.3]);
ylim([-1.3 1.3]);
grid on;
axis square;
ax = gca;
ax.FontSize = 16;
xlabel('\boldmath$\hat{s}_1$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('\boldmath$\hat{s}_2$', 'Interpreter', 'latex', 'FontSize', 20, 'Rotation', 0);

vertices = round(P_x_set_tight.V, 2);
figure();
fill(vertices(:,1), vertices(:,2), [0 0.4470 0.7410], 'FaceAlpha', 0.4);
xlim([-1.3 1.3]);
ylim([-1.3 1.3]);
grid on;
axis square;
text(vertices(1,1), vertices(1,2), ['(' num2str(vertices(1,1)) ',' num2str(vertices(1,2)) ')'], 'HorizontalAlignment','left', 'FontSize', 16);
text(vertices(2,1), vertices(2,2), ['(' num2str(vertices(2,1)) ',' num2str(vertices(2,2)) ')'], 'HorizontalAlignment','right', 'FontSize', 16);
text(vertices(3,1), vertices(3,2), ['(' num2str(vertices(3,1)) ',' num2str(vertices(3,2)) ')'], 'HorizontalAlignment','right', 'FontSize', 16);
text(vertices(4,1), vertices(4,2), ['(' num2str(vertices(4,1)) ',' num2str(vertices(4,2)) ')'], 'HorizontalAlignment','left', 'FontSize', 16);
ax = gca;
ax.FontSize = 16;
xlabel('\boldmath$\bar{x}_1$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('\boldmath$\bar{x}_2$', 'Interpreter', 'latex', 'FontSize', 20, 'Rotation', 0);

vertices = round(P_rect_incr.V, 2);
figure();
fill(vertices(:,1), vertices(:,2), [0 0.4470 0.7410], 'FaceAlpha', 0.4);
xlim([-1.3 1.3]);
ylim([-1.3 1.3]);
grid on;
axis square;
text(vertices(1,1), vertices(1,2), ['(' num2str(vertices(1,1)) ',' num2str(vertices(1,2)) ')'], 'HorizontalAlignment','right', 'FontSize', 16);
text(vertices(2,1), vertices(2,2), ['(' num2str(vertices(2,1)) ',' num2str(vertices(2,2)) ')'], 'HorizontalAlignment','left', 'FontSize', 16);
text(vertices(3,1), vertices(3,2), ['(' num2str(vertices(3,1)) ',' num2str(vertices(3,2)) ')'], 'HorizontalAlignment','left', 'FontSize', 16);
text(vertices(4,1), vertices(4,2), ['(' num2str(vertices(4,1)) ',' num2str(vertices(4,2)) ')'], 'HorizontalAlignment','right', 'FontSize', 16);
ax = gca;
ax.FontSize = 16;
xlabel('\boldmath$\hat{s}_1$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('\boldmath$\hat{s}_2$', 'Interpreter', 'latex', 'FontSize', 18, 'Rotation', 0);

vertices = round(P_x_set_tight_rough.V, 2);
figure();
fill(vertices(:,1), vertices(:,2), [0 0.4470 0.7410], 'FaceAlpha', 0.4);
xlim([-1.3 1.3]);
ylim([-1.3 1.3]);
grid on;
axis square;
text(vertices(1,1), vertices(1,2), ['(' num2str(vertices(1,1)) ',' num2str(vertices(1,2)) ')'], 'HorizontalAlignment','left', 'FontSize', 16);
text(vertices(2,1), vertices(2,2), ['(' num2str(vertices(2,1)) ',' num2str(vertices(2,2)) ')'], 'HorizontalAlignment','right', 'FontSize', 16);
text(vertices(3,1), vertices(3,2), ['(' num2str(vertices(3,1)) ',' num2str(vertices(3,2)) ')'], 'HorizontalAlignment','right', 'FontSize', 16);
text(vertices(4,1), vertices(4,2), ['(' num2str(vertices(4,1)) ',' num2str(vertices(4,2)) ')'], 'HorizontalAlignment','left', 'FontSize', 16);
ax = gca;
ax.FontSize = 16;
xlabel('\boldmath$x_1$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('\boldmath$x_2$', 'Interpreter', 'latex', 'FontSize', 18, 'Rotation', 0);

% Define \bar{X} in terms of Ax <= b
A_c = P_x_set_tight.A;
b_c = P_x_set_tight.b;



% Define original constraint set U
% |u| <= 1
u_max = 1;
u_min = -1;

% Compute upper and lower bound of KS(infty)
sol_max = linprog(-K.',A_s,b_s);
max_ks = K*sol_max;

sol_min = linprog(K.',A_s,b_s);
min_ks = K*sol_min;

% Compute tightend constraint set \bar{U}
% \bar{U} = U - KS(infty)
u_tight_max = u_max - max_ks;
u_tight_min = u_min - min_ks;


% Store tightened constraints for use outside this script
save('tightened_constraints.mat', 'A_c', 'b_c', 'u_tight_max', 'u_tight_min');


% fh = findall(0,'Type','Figure');
% set(findall(fh, '-property', 'fontsize'), 'fontsize', 16);



function x = simulate_sys(A, B, u, x_0, n)
    x = [x_0 zeros(n, size(u,2))];
    for i = 2:1+size(u,2)
        x(:, i) = A*x(:,i-1) + B*u(:,i-1);
    end
end

function [x, u] = simulate_uncertain_sys(A, B, u_n, x_n, x_0, n, m, K, w)
    x = [x_0 zeros(n, size(u_n,2))];
    u = zeros(m, size(u_n,2));
    for i = 2:1+size(u_n,2)
        u(:,i-1) = u_n(:,i-1) + K*(x(:,i-1)-x_n(:,i-1));
        x(:,i) = A*x(:,i-1) + B*u(:,i-1) + w(:,i-1);
    end
end



