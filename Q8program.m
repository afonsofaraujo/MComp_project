%Q8
close all;
clear;
clc;

%[coordout, connectivityData] = readNXData('elements Q8.txt', 'nos Q8.txt');
[coordout, connectivityData] = readNXData('el8.txt', 'nd8.txt');
disp('Data loaded...');
[fronteira, B1, B2, B3, B4] = identifyBoundary(coordout, 0.3);
disp('Boundary nodes identified...');
disp('Assembly...');
[Kg, fg] = assembleGlobalMatrixAndForce(coordout, connectivityData);
disp('Boundary conditions...');
U = 2.5;
[Kg, fg] = applyBoundaryConditions(Kg, fg, B1, B2, B4, coordout, U);
disp('Boundary conditions applied...');
disp('Solution...');
u = solveSystem(Kg, fg);

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
