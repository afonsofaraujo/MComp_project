%Solução para elementos de 4 nós - Q4
close all;
clear;
clc;

%Leitura do ficheiro msh

[coordout, connectivityData] = readNXData('Elementscomsimetria.txt', 'Nodescomsimetria.txt');
disp('Data loaded...');

%Definição da fronteira

[fronteira, B1, B2, B3, B4] = identifyBoundary(coordout, 0.85);
disp('Boundary nodes identified...');
disp('Assembly...');

[Kg, fg] = assembleGlobalMatrixAndForce(coordout, connectivityData); 
disp('Global matrix assembled...');
disp('Boundary conditions...');
U = 2.5;    % velocidade de entrada em m/s
[Kg, fg] = applyBoundaryConditions(Kg, fg, B1, B2, B4, coordout, U);
disp('Boundary conditions applied...');
disp('Solution...');
u = solveSystem(Kg, fg);
disp('System solved...');


%%%%%%%%%%%%%%
disp('Post-processing...');

elementType = 'QUAD4';
Nels = size(connectivityData, 1);
Nnds = length(coordout(:, 2));
coordx = coordout(:,2);
coordy = coordout(:,3);

[xcentroid, ycentroid] = computeCentroids(connectivityData, coordx, coordy, elementType);
pressure = calculatePressure(connectivityData, coordx, coordy, u, elementType);
%[fineX, fineY, fineU, fineP] = interpolateAndMask(coordx, coordy, u, pressure, xcentroid, ycentroid, fronteira, 2);
[fineXY, fineU, fineP] = interpolateOption(coordx, coordy, u, xcentroid, ycentroid, pressure, fronteira, 1);
[vx, vy] = calculateCentroidsVelocity(connectivityData, coordx, coordy, u, elementType);
[Res, xint, yint, vxint, vyint] = calculateVelocityAtIntegrationPoints(connectivityData, coordx, coordy, u, 4, elementType);



surf(fineXY(:,1),fineXY(:,2),fineU(:))



% Scatter plot for fineU
figure;
scatter(fineXY(:, 1), fineXY(:, 2), 50, fineU, 'filled');
colorbar;
title('Fine U Distribution');
xlabel('X-axis');
ylabel('Y-axis');
axis equal;

% Scatter plot for fineP
figure;
scatter(fineXY(:, 1), fineXY(:, 2), 50, fineP, 'filled');
colorbar;
title('Fine P Distribution');
xlabel('X-axis');
ylabel('Y-axis');
axis equal;














%%%%%%%%%%%%%% streamlines

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

%%%%%%%%%%%%%% pressure

figure;
contour(fineX, fineY, fineP, 10,'LineColor','k',"ShowText",true,"LabelFormat","%0.2f bar",'LabelSpacing', 400);
hold on;
plot(fronteira(:,1), fronteira(:,2), 'Color', 'k', 'LineWidth', 1);    % Plot boundary
title('P');
xlabel('X-axis');
ylabel('Y-axis');
axis equal;
hold off;

%%%%%%%%%%%%%% velocity

figure;

subplot(2, 1, 1);
quiver(xint, yint, vxint, vyint);
axis equal;
title('Pontos de integração');
xlabel('X-axis');
ylabel('Y-axis');
hold off;

subplot(2, 1, 2);
quiver(xcentroid, ycentroid, vx, vy);
axis equal;
title('Centróides');
xlabel('X-axis');
ylabel('Y-axis');
hold off;

%%%%%%%%%%%%%% MIRROR pressure

% Create a mirrored version of pressure along Y = 0.9 axis
mirroredFineP = flipud(fineP);

% Create meshgrid for mirrored coordinates
fineXmir = linspace(min(coordx), max(coordx), 2 * (length(coordx) - 1) + 1);
fineYmir = linspace(max(coordy), 2 * max(coordy), 2 * (length(coordy) - 1) + 1);
[fineXmir, fineYmir] = meshgrid(fineXmir, fineYmir);

% Combine original and mirrored data
combinedFineX = [fineX, fineXmir];
combinedFineY = [fineY, fineYmir];
combinedFineP = [fineP, mirroredFineP];

fill3(combinedFineX, combinedFineY,combinedFineP,combinedFineP);
% Plot solution
figure
for i=1:Nels
    edofs=[connectivityData(i,:)];
    fill3(coordx(edofs),coordy(edofs),u(edofs),u(edofs));hold on
end
plot(coordx,coordy,'ro');
hold off;




% Plot contours for the entire figure
figure;

contour(combinedFineX, combinedFineY, combinedFineP, 10, 'LineColor', 'k', 'ShowText', true, 'LabelFormat', '%0.2f bar', 'LabelSpacing', 400);
hold on;

% Plot boundary
B5 = [B1(:,1), 2*max(coordy)-B1(:,2)];
B6 = [B3(:,1), 2*max(coordy)-B3(:,2)];
B7 = [B2(:,1), 2*max(coordy)-B2(:,2)];
plot(B1(:,1), B1(:,2),'k', 'LineWidth', 1),
plot(B2(:,1), B2(:,2),'k', 'LineWidth', 1);
plot(B3(:,1), B3(:,2),'k', 'LineWidth', 1);
plot(B5(:,1), B5(:,2),'k', 'LineWidth', 1);
plot(B6(:,1), B6(:,2),'k', 'LineWidth', 1);
plot(B7(:,1), B7(:,2),'k', 'LineWidth', 1);

title('P');
xlabel('X-axis');
ylabel('Y-axis');
axis equal;
legend('Combined Boundary');
hold off;







%%%%%%%%% Mirror stream lines

figure;
contourf(fineX, fineY, fineU, 20);
hold on;
plot(fronteira(:,1), fronteira(:,2), 'Color', 'k', 'LineWidth', 1);    % Plot boundary
title('U');
xlabel('X-axis');
ylabel('Y-axis');
axis equal;
hold off;





























% % Original pressure plot
% figure;
% 
% contour(fineX, fineY, fineP, 10,'LineColor','k',"ShowText",true,"LabelFormat","%0.2f bar",'LabelSpacing', 400);
% hold on;
% plot(fronteira(:,1), fronteira(:,2), 'Color', 'k', 'LineWidth', 1);    % Plot boundary
% title('Original Pressure');
% xlabel('X-axis');
% ylabel('Y-axis');
% axis equal;
% hold off;

% % Create a mirrored version of pressure along y = 3.8 axis
% mirroredPressure = flipud(fineP);
% 
% % Plot mirrored pressure along y = 3.8 axis
% figure;
% contour(fineX, fineY, mirroredPressure, 10,'LineColor','k',"ShowText",true,"LabelFormat","%0.2f bar",'LabelSpacing', 400);
% hold on;
% plot(fronteira(:,1), fronteira(:,2), 'Color', 'k', 'LineWidth', 1);    % Plot boundary
% title('Mirrored Pressure along y = 3.8 axis');
% xlabel('X-axis');
% ylabel('Y-axis');
% axis equal;
% hold off;

