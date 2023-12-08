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

half = 0; % For half piece half = 1, whole half = 0,

if half

    [xcentroid, ycentroid] = computeCentroids(connectivityData, coordx, coordy, elementType);
    pressure = calculatePressure(connectivityData, coordx, coordy, u, elementType);
    [fineX, fineY, fineU, fineP] = interpolateAndMask(coordx, coordy, u, pressure, xcentroid, ycentroid, fronteira, 2);
    [vx, vy] = calculateCentroidsVelocity(connectivityData, coordx, coordy, u, elementType);
    [Res, xint, yint, vxint, vyint] = calculateVelocityAtIntegrationPoints(connectivityData, coordx, coordy, u, 4, elementType);
    
    % Streamlines
    figure, subplot(2, 1, 1);
    contourf(fineX, fineY, fineU, 20), hold on;
    plot(fronteira(:,1), fronteira(:,2), 'Color', 'k', 'LineWidth', 1);
    title('U'), xlabel('X-axis'), ylabel('Y-axis'), axis equal;hold off;
    
    subplot(2, 1, 2);
    contour(fineX, fineY, fineU, 20,'LineColor','k'), hold on;
    plot(fronteira(:,1), fronteira(:,2), 'Color', 'k', 'LineWidth', 1);
    title('U'), xlabel('X-axis'), ylabel('Y-axis'), axis equal, hold off;
    
    % Pressure
    figure;
    contour(fineX, fineY, fineP, 10,'LineColor','k',"ShowText",true,"LabelFormat","%0.2f bar",'LabelSpacing', 400), hold on;
    plot(fronteira(:,1), fronteira(:,2), 'Color', 'k', 'LineWidth', 1);
    title('P'), xlabel('X-axis'), ylabel('Y-axis'), axis equal, hold off;
    
    % Velocity
    figure, subplot(2, 1, 1);
    quiver(xint, yint, vxint, vyint);
    title('Pontos de integração'), xlabel('X-axis'), ylabel('Y-axis'), axis equal, hold off;

    subplot(2, 1, 2);
    quiver(xcentroid, ycentroid, vx, vy);
    title('Centróides'), xlabel('X-axis'), ylabel('Y-axis'), axis equal, hold off;
else
    pressure = calculatePressure(connectivityData, coordx, coordy, u, elementType);
    px = [coordx(:); coordx(:)];
    py = [coordy(:); 2*max(coordy(:))-coordy(:)];
    [pxpy, uniqueIndices] = unique([px, py], 'rows');
    ut = [u; u];
    ut = ut(uniqueIndices, :);  % Keep only the rows corresponding to unique points
    px = pxpy(:,1);
    py = pxpy(:,2);
    B5 = flipud([B3(:,1), 2*max(coordy)-B3(:,2)]);
    B6 = flipud([B2(:,1), 2*max(coordy)-B2(:,2)]);
    B7 = [B1(:,1), 2*max(coordy)-B1(:,2)];
    fronteirat = [B1;B2;B3;B5;B6;B7];
    [xcentroid, ycentroid] = computeCentroids(connectivityData, coordx, coordy, elementType);
    pxcentroid = [xcentroid; xcentroid];
    pycentroid = [ycentroid; 2*max(coordy)-ycentroid];
    pressuret = [pressure; pressure];
    
    [fineX, fineY, fineUt, finePt] = interpolateAndMask(px, py, ut, pressuret, pxcentroid, pycentroid, fronteirat, 2);
    
    % Stream lines
    figure, subplot(2, 1, 1);
    contourf(fineX, fineY, fineUt, 20), hold on;
    plot([min(px), max(px)], [0.9, 0.9], 'k-', 'LineWidth', 0.5);
    plot(fronteirat(:,1), fronteirat(:,2), 'Color', 'k', 'LineWidth', 1);
    title('U'), xlabel('X-axis'), ylabel('Y-axis'), axis equal, axis off, hold off;
    subplot(2, 1, 2);
    contour(fineX, fineY, fineUt, 20,'LineColor','k'), hold on;
    plot([min(px), max(px)], [0.9, 0.9], 'k-', 'LineWidth', 0.5);
    plot(fronteirat(:,1), fronteirat(:,2), 'Color', 'k', 'LineWidth', 1);
    title('U'), xlabel('X-axis'), ylabel('Y-axis'), axis equal, axis off, hold off;
    
    % Pressure
    figure;
    contour(fineX, fineY, finePt, 10,'LineColor','k',"ShowText",true,"LabelFormat","%0.2f bar",'LabelSpacing', 400), hold on;
    plot(fronteirat(:,1), fronteirat(:,2), 'Color', 'k', 'LineWidth', 1);
    title('P'), xlabel('X-axis'), ylabel('Y-axis'), axis equal, hold off;
end