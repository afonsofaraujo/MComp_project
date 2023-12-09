% Solução para elementos de 8 nós - Q8
close all;
clear, clc;

elementType = 'QUAD8';
boundaryParameter = 0.3;
elementsFileName = 'el8.txt';
nodesFileName = 'nd8.txt';
U = 2.5;    % velocidade de entrada em m/s

% Leitura do ficheiro
[coordout, connectivityData] = readNXData(elementsFileName, nodesFileName);
disp('Data loaded...');

% Definição da fronteira
[fronteira, B1, B2, B3, B4] = identifyBoundary(coordout, boundaryParameter);
disp('Boundary nodes identified...');

% Assemblagem
[Kg, fg] = assembleGlobalMatrixAndForce(coordout, connectivityData); 
disp('Global matrix assembled...');

% Condições de fronteira
[Kg, fg] = applyBoundaryConditions(Kg, fg, B1, B2, B4, coordout, U);
disp('Boundary conditions applied...');

% Resolução do sistema
u = solveSystem(Kg, fg);
disp('System solved...');

% Pós-processamento
disp('Post-processing...');

Nels = size(connectivityData, 1);
Nnds = length(coordout(:, 2));
coordx = coordout(:,2);
coordy = coordout(:,3);

half = 1; % For half piece half = 1, whole half = 0,

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
    title('Pontos de integração'), xlabel('X-axis'), ylabel('Y-axis'), axis equal;

    subplot(2, 1, 2);
    quiver(xcentroid, ycentroid, vx, vy);
    title('Centróides'), xlabel('X-axis'), ylabel('Y-axis'), axis equal;
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
    [Res, xint, yint, vxint, vyint] = calculateVelocityAtIntegrationPoints(connectivityData, coordx, coordy, u, 4, elementType);
    xintt = [xint; xint];
    yintt = [yint; 2*max(coordy)-yint];
    vxintt = [vxint; vxint];
    vyintt = [vyint; -vyint];

    [vx, vy] = calculateCentroidsVelocity(connectivityData, coordx, coordy, u, elementType);
    vxt = [vx; vx];
    vyt = [vy; -vy];
    
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

    % Velocity
    figure, subplot(2, 1, 1);
    quiver(xintt, yintt, vxintt, vyintt), hold on;
    B6 = [B2(:,1), 2*max(coordy)-B2(:,2)];
    plot(B2(:,1), B2(:,2),'k', 'LineWidth', 0.5);
    plot(B6(:,1), B6(:,2),'k', 'LineWidth', 0.5);
    %plot(fronteirat(:,1), fronteirat(:,2), 'Color', 'k', 'LineWidth', 0.5);
    title('Pontos de integração'), xlabel('X-axis'), ylabel('Y-axis'), axis equal, hold off;
    subplot(2, 1, 2);
    quiver(pxcentroid, pycentroid, vxt, vyt), hold on;
    plot(B2(:,1), B2(:,2),'k', 'LineWidth', 0.5);
    plot(B6(:,1), B6(:,2),'k', 'LineWidth', 0.5);
    %plot(fronteirat(:,1), fronteirat(:,2), 'Color', 'k', 'LineWidth', 0.5);
    title('Centróides'), xlabel('X-axis'), ylabel('Y-axis'), axis equal, hold off;
end
%%%%%%%
disp('Post-processing...');

elementType = 'QUAD8';
Nels = size(connectivityData, 1);
Nnds = length(coordout(:, 2));
coordx = coordout(:,2);
coordy = coordout(:,3);

[xcentroid, ycentroid] = computeCentroids(connectivityData, coordx, coordy, elementType);
pressure = calculatePressure(connectivityData, coordx, coordy, u, elementType);
[fineX, fineY, fineU, fineP] = interpolateAndMask(coordx, coordy, u, pressure, xcentroid, ycentroid, fronteira, 1); % mesh density can be change (last parameter)
[vx, vy] = calculateCentroidsVelocity(connectivityData, coordx, coordy, u, elementType);
[Res, xint, yint, vxint, vyint] = calculateVelocityAtIntegrationPoints(connectivityData, coordx, coordy, u, 4, elementType);


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

%%%%%%%%%%%%%

% Plot elements
connectivityPatch = zeros(size(connectivityData));
for i = 1:length(connectivityData)
    connectivityPatch(i,:) = [connectivityData(i,1), connectivityData(i,5),...
    connectivityData(i,2), connectivityData(i,6), connectivityData(i,3),...
    connectivityData(i,7), connectivityData(i,4), connectivityData(i,8)];
end

figure;
hold on;
for i = 1:Nels
    currentConnectivity = connectivityPatch(i, :);
    x = coordx(currentConnectivity);
    y = coordy(currentConnectivity);
    patch(x, y, 'b', 'FaceAlpha', 0.5);
    for j = 1:length(currentConnectivity)
        scatter(x(j), y(j), 20, 'r', 'filled');
        text(x(j), y(j), num2str(currentConnectivity(j)), 'Color', 'k', 'FontSize', 8);
    end
end
axis equal;
axis off;
hold off;






% % Plot boundary
% plot(B1(:,1), B1(:,2)), hold on;
% plot(B2(:,1), B2(:,2));
% plot(B3(:,1), B3(:,2));
% plot(B4(:,1), B4(:,2)); 
% legend(["B1","B2","B3","B4"]);
% hold off;
% plot(fronteira(:,1), fronteira(:,2))

% % Plot solution
% figure
% for i=1:Nels
%     edofs=[connectivityPatch(i,:)];
%     fill3 (coordx(edofs),coordy(edofs),u(edofs),u(edofs));hold on
%     plot(coordx(edofs),coordy(edofs),'r');hold on
% end
% plot(coordx,coordy,'ro');
% hold off;
