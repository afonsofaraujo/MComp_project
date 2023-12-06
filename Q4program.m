%Q4
close all;
clear;
clc;

[coordout, connectivityData] = readNXData('Elementscomsimetria.txt', 'Nodescomsimetria.txt');
disp('Data loaded...');
[fronteira, B1, B2, B3, B4] = identifyBoundary(coordout);
disp('Assembly...');
[Kg, fg] = assembleGlobalMatrixAndForce(coordout, connectivityData); 
disp('Global matrix assembled...');
disp('Boundary conditions...');
U = 2.5;
[Kg, fg] = applyBoundaryConditions(Kg, fg, B1, B2, B4, coordout, U);
disp('Solution...');
u = solveSystem(Kg, fg);
disp('Post-processing...');

Nels = size(connectivityData, 1);
Nnds = length(coordout(:, 2));
coordx = coordout(:,2);
coordy = coordout(:,3);

constant = 104450;  % 101325+0.5*2.5^2*1000
rho = 1000;         % 1000 water

pressure = zeros(Nels, 1);
fluxArray = zeros(Nels, 2); % Preallocate arrays for storing calculated values
centroidArray = zeros(Nels, 2);
for i = 1:Nels
    edofs = connectivityData(i, :);
    XN(1:4, 1) = coordx(edofs); % Extract coordinates
    XN(1:4, 2) = coordy(edofs);
    csi = 0;
    eta = 0;
    [B, psi, Detj] = Shape_N_Der4(XN, csi, eta);    % Calculate shape functions and derivatives
    uint = psi' * u(edofs);
    xpint = XN' * psi;    % Position (x, y) of the centroid
    gradu = B' * u(edofs);
    fluxu = -gradu;
    fluxArray(i, :) = fluxu;        % store flux values
    centroidArray(i, :) = xpint;    % centroid coordinates
    pressure(i) = (constant - 0.5*rho*(sqrt(fluxu(1)^2+fluxu(2)^2))^2)*10^(-5);  % pressure in bar
end
 
%%%%%%%%%%%%%%%%

% Define finer mesh grid
fineMeshDensity = 2;  % Adjust this value based on your needs
fineX = linspace(min(coordx), max(coordx), fineMeshDensity * (length(coordx)-1) + 1);
fineY = linspace(min(coordy), max(coordy), fineMeshDensity * (length(coordy)-1) + 1);
[fineX, fineY] = meshgrid(fineX, fineY);

fineU = griddata(coordx, coordy, u, fineX, fineY);  % Interpolate nodal displacements on the finer mesh
fineP = griddata(centroidArray(:,1), centroidArray(:,2), pressure, fineX, fineY);
inBoundary = inpolygon(fineX, fineY, fronteira(:,1), fronteira(:,2));
fineU(~inBoundary) = NaN;   % Set values outside the boundary to NaN
fineP(~inBoundary) = NaN;   % Set values outside the boundary to NaN

disp(['Min fineU: ', num2str(min(fineU(:)))]);
disp(['Max fineU: ', num2str(max(fineU(:)))]);

figure;
subplot(2, 1, 1);
contourf(fineX, fineY, fineU, 20);
hold on;
plot(fronteira(:,1), fronteira(:,2), 'Color', 'k', 'LineWidth', 1);    % Plot boundary
title('U');
xlabel('X-axis');
ylabel('Y-axis');
axis equal;
hold off;

subplot(2, 1, 2);
contour(fineX, fineY, fineU, 20,'LineColor','k');
hold on;
plot(fronteira(:,1), fronteira(:,2), 'Color', 'k', 'LineWidth', 1);    % Plot boundary
title('U');
xlabel('X-axis');
ylabel('Y-axis');
axis equal;
hold off;

%%%%%%%%%%%%%%

figure;
contour(fineX, fineY, fineP, 10,'LineColor','k',"ShowText",true,"LabelFormat","%0.2f bar",'LabelSpacing', 400);
hold on;
plot(fronteira(:,1), fronteira(:,2), 'Color', 'k', 'LineWidth', 1);    % Plot boundary
title('P');
xlabel('X-axis');
ylabel('Y-axis');
axis equal;
hold off;

%%%%%%%%%%%%%%%%%%

Res = 0;
arrow_x = [];   % Initialize arrays to store results
arrow_y = [];
arrow_u = [];
arrow_v = [];
for i = 1:Nels
    edofs = connectivityData(i,:);
    XN(1:4,1) = coordx(edofs);
    XN(1:4,2) = coordy(edofs);
    csi = 0;
    eta = 0;
    nip = 4;
    [xp, wp] = Genip2DQ(nip);
    for ip = 1:nip
        csi = xp(ip,1);
        eta = xp(ip,2);
        [B, psi, Detj] = Shape_N_Der4(XN, csi, eta);
        uint = psi' * u(edofs);
        Res = Res + uint * Detj * wp(ip);
        xpint = XN' * psi;
        gradu = B' * u(edofs);
        fluxu = -gradu;
        arrow_x = [arrow_x; xpint(1)];  % Store results for quiver plot
        arrow_y = [arrow_y; xpint(2)];
        arrow_u = [arrow_u; fluxu(1)];
        arrow_v = [arrow_v; fluxu(2)];
    end
end
rotationMatrix = [cosd(90) -sind(90); sind(90) cosd(90)];
rotatedVectors = rotationMatrix * [arrow_u(:)'; arrow_v(:)'];   % Apply the rotation to the vectors
u2 = reshape(rotatedVectors(1, :), size(arrow_x));    % Reshape the rotated vectors back to the original grid
v2 = reshape(rotatedVectors(2, :), size(arrow_y));

%%%%%%%%%%%%

fluxArray = zeros(Nels, 2); % Preallocate arrays for storing calculated values
centroidArray = zeros(Nels, 2);
for i = 1:Nels
    edofs = connectivityData(i, :);
    XN(1:4, 1) = coordx(edofs); % Extract coordinates
    XN(1:4, 2) = coordy(edofs);
    csi = 0;
    eta = 0;
    [B, psi, Detj] = Shape_N_Der4(XN, csi, eta);    % Calculate shape functions and derivatives
    uint = psi' * u(edofs);
    xpint = XN' * psi;    % Position (x, y) of the centroid
    gradu = B' * u(edofs);
    fluxu = -gradu;
    fluxArray(i, :) = fluxu;        % store flux values
    centroidArray(i, :) = xpint;    % centroid coordinates
end
rotationMatrix = [cosd(90) -sind(90); sind(90) cosd(90)];
rotatedVectors = rotationMatrix * [fluxArray(:, 1)'; fluxArray(:, 2)'];   % Apply the rotation to the vectors
vec_x = rotatedVectors(1, :)';    % Reshape the rotated vectors back to the original grid
vec_y = rotatedVectors(2, :)';
% for i = 1:Nels
%     currentConnectivity = connectivityData(i, :);
%     x = coordx(currentConnectivity);
%     y = coordy(currentConnectivity);
%     patch(x, y, norm(fluxArray(i, :))),hold on;
% end

%%%%%%%%%%%%%%%%%%%%%

figure;

subplot(2, 1, 1);
quiver(arrow_x, arrow_y, u2, v2);
axis equal;
title('Pontos de integração');
xlabel('X-axis');
ylabel('Y-axis');
hold off;

subplot(2, 1, 2);
quiver(centroidArray(:, 1), centroidArray(:, 2), vec_x, vec_y);       % plot velocity arrows
axis equal;
title('Centróides');
xlabel('X-axis');
ylabel('Y-axis');
hold off;
