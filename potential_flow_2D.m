% See also http://www.particleincell.com/blog/2012/matlab-fem/
% 2D Potential flow
close all; clear; clc
 
% Mesh loading, set parameters
load mesh_vessel.mat
number_of_nodes = size(p, 2);
number_of_elements = size(t, 1);
 
% Initialization of K and F
K = zeros(number_of_nodes);
F = zeros(number_of_nodes, 1);
 
% Calculation of Ke & assembly of K 
for element = 1 : number_of_elements
    nodes = t(element, :);
    P = [ones(1, 3); p(:, nodes)];
    C = inv(P);
    area_of_element = abs(det(P))/2;
    grads_phis = C(:, 2:3);
    xy_mean = mean(p(:, nodes), 2);
    Ke = grads_phis * grads_phis' * area_of_element;
    K(nodes, nodes) = K(nodes, nodes) + Ke;
end

% Neumann boundary
q = -1; % initial velocity
t_Neumann = [1 14; 14 2; 7 25; 25 6];
for element = 1 : size(t_Neumann, 1)
    nodes = t_Neumann(element, :);
    P = p(:, nodes);
    length_of_element = norm(diff(P, 1, 2));
    xy_mean = mean(p(:, nodes), 2);
    Fe = 1/2 * q * length_of_element * [1; 1];
    F(nodes) = F(nodes) + Fe;
end
 
% Dirichlet boundary
Dirichlet = [4 24 5];
K(Dirichlet, :) = 0;
K(Dirichlet, Dirichlet) = eye(numel(Dirichlet));
F(Dirichlet) = 0;
 
% Solve
phi = K \ F;
 
% Plotting
trisurf(t, p(1, :), p(2, :), phi, 'edgecolor', ...
    0.5 * [1 1 1], 'facecolor', 'none');
 
% Calculation of velocity vectors at center of elements
Z = [];
V = [];
magnitudes_nodal_velocities = zeros(number_of_nodes, 1);
for element = 1 : number_of_elements
    nodes = t(element, :);
    P = [ones(1, 3); p(:, nodes)];
    C = inv(P);
    grads_phis = C(:, 2:3);
    xy_mean = mean(p(:, nodes), 2);
    Z = [Z; xy_mean'];
    phi_at_nodes = phi(nodes,:);
    v_e = grads_phis' * phi_at_nodes;
    V = [V; v_e'];
    magnitudes_nodal_velocities(nodes) = magnitudes_nodal_velocities(nodes) + norm(v_e);
end

% Determine nodal velocity magnitudes
node_frequencies = [];
for i = 1 : number_of_nodes
    text(p(1, i), p(2, i), int2str(i), 'FontName', 'Consolas', ...
        'FontSize', 10, 'HorizontalAlignment', 'Center', ...
        'FontWeight', 'bold', 'Color', 0.75 * [1 1 1])
    node_frequencies = [node_frequencies; numel(find(t == i))];
end
magnitudes_nodal_velocities = magnitudes_nodal_velocities ./ node_frequencies;

% Plot velocity field
set(gcf, 'Color', 'white')
view(2)
hold on
trisurf(t, p(1, :), p(2, :), 0 * p(1, :), magnitudes_nodal_velocities, 'edgecolor', 'w', 'facecolor', 'interp');
quiver(Z(:, 1), Z(:, 2), V(:, 1), V(:, 2), 'y', 'LineWidth', 1.5)
colormap jet
caxis([min(magnitudes_nodal_velocities), max(magnitudes_nodal_velocities)])
axis equal
axis off
 
% Calculate streamlines
x = Z(:, 1); 
y = Z(:, 2); 
u = V(:, 1);
v = V(:, 2);

% Interpolate
U = scatteredInterpolant(x, y, u);
V = scatteredInterpolant(x, y, v);
x_interp = linspace(min(x), max(x), 50);
y_interp = linspace(min(y), max(y), 50);
[X, Y] = meshgrid(x_interp, y_interp);
u_interp = U(X, Y);
v_interp = V(X, Y);

% Startpositions particles
p_1 = mean([p(:, 1), p(:, 14)], 2); p_1(1) = p_1(1) + 1e-1; p_1(2) = p_1(2) - 1e-1;
p_2 = mean([p(:, 14), p(:, 2)], 2); p_2(1) = p_2(1) + 1e-1; p_2(2) = p_2(2) - 1e-1;
p_3 = mean([p(:, 6), p(:, 25)], 2); p_3(1) = p_3(1) + 1e-1;
p_4 = mean([p(:, 25), p(:, 7)], 2); p_4(1) = p_4(1) + 1e-1;
x_ = [p_1(1) p_2(1) p_3(1) p_4(1)];
y_ = [p_1(2) p_2(2) p_3(2) p_4(2)];

% Streamline generation
Iverts = [];
for i = 1 : numel(x_)
    h(i) = streamline(stream2(X, Y, u_interp, v_interp, x_(i), y_(i)));
    set(h, 'Color', [1, 0, 0], 'LineWidth', 1.5)   
    verts = stream2(X, Y, u_interp, v_interp, x_(i), y_(i));
    iverts = interpstreamspeed(X, Y, u_interp, v_interp, verts, 0.002);
    Iverts = [Iverts, iverts];
end

% Animation
streamparticles(Iverts, 20, ...
    'Animate', 20, ...
    'ParticleAlignment', 'on', ...
    'Marker', 'o', ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'r', ...
    'MarkerSize', 4.5);