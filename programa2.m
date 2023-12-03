%PROGRAMA2
close all; clear; clc
%[coordout, connectivityData] = readNXData('Elementscasosimples.txt', 'Nodescasosimples.txt');
[coordout, connectivityData] = readNXData('Elementscomsimetria.txt', 'Nodescomsimetria.txt');

disp('Data loaded...');
coordx = coordout(:,2);
coordy = coordout(:,3);
k = boundary(coordx,coordy,0.85);
fronteira = [coordx(k),coordy(k)];
ymax = max(fronteira(:,2));
xmax = max(fronteira(:,1));
xmin = min(fronteira(:,1));
B1 = [];    % esquerda
B2 = [];    % baixo
B3 = [];    % direita
B4 = [];    % cima
s = [];     % índices já utilizados
for i = 1:size(fronteira,1)
    if fronteira(i,2) == ymax
        B4 = [B4; fronteira(i,:)];
        s = [s;i];
    end
    if fronteira(i,1) == xmin
        B1 = [B1; fronteira(i,:)];
        s = [s;i];
    end
    if fronteira(i,1) == xmax
        B3 = [B3; fronteira(i,:)];
        s = [s;i];
    end
end
s = unique(s);  % retirar índices repetidos (cantos)
for i = 1:size(fronteira,1)
    if ~ismember(i,s)
        B2 = [B2; fronteira(i,:)]; % B2 contém os pontos da fronteira não utilizados
    end
end
B2 = [B2; xmax min(B3(:,2))];   % adicionar ponto final 
B2 = [0 0; B2];                 % adicionar ponto inicial

disp('Assembly...');
Nels = size(connectivityData, 1);
Nnds = length(coordx);
Kg = zeros(Nnds,Nnds); % inicialização
fg = zeros(Nnds,1);
for i=1:Nels    % Carregamento no elemento
    edofs=[connectivityData(i, :)];
    XN(1:4,1)=coordx(edofs);
    XN(1:4,2)=coordy(edofs);
    [Ke, fe]=Elem_Quad4(XN,0);  % carregamento 0
    Kg(edofs,edofs)= Kg(edofs,edofs) + Ke;
    fg(edofs,1)= fg(edofs,1) + fe;
end

disp('Boundary conditions...');
boom = 1.0e+14;
U = 2.5*1000;   % velocidade de entrada em mm/s (porque a malha está em mm)
for i = 1:size(B1,1)    % B1
    j = find(coordx == B1(i, 1) & coordy == B1(i, 2));
    Kg(j, j) = boom;
    fg(j) = boom*coordy(j)*U; % stream function = yU
end
for i = 1:size(B2,1)    % B2
    j = find(coordx == B2(i, 1) & coordy == B2(i, 2)); % Find the corresponding node
    Kg(j, j) = boom;
    fg(j) = 0; % stream function = 0
end
for i = 1:size(B4,1)    % B4
    j = find(coordx == B4(i, 1) & coordy == B4(i, 2));
    Kg(j, j) = boom;
    fg(j) = boom*ymax*U; % stream function = y_maxU
end

disp('Solution...');
u=Kg\fg;    % Resolução do sistema
disp('Post-processing...');

constant = 10^12;
rho = 1;

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
    pressure(i) = constant - 0.5*rho*(sqrt(fluxu(1)^2+fluxu(2)^2)*1000)^2;  % P + 0.5*ρ*v^2 = constant
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

%%%%%%%%%%%%%%%

figure;

contour(fineX, fineY, fineP, 10,'LineColor','k');
%contourf(fineX, fineY, fineP, 10);
colorbar;
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

%% Plot elements
% figure;
% hold on;
% for i = 1:Nels
%     currentConnectivity = connectivityData(i, :);
%     x = coordx(currentConnectivity);
%     y = coordy(currentConnectivity);
%     patch(x, y, 'b', 'FaceAlpha', 0.5);
% end
% hold off;
% axis equal;
% grid on;
% xlabel('X-axis');
% ylabel('Y-axis');
% title('Quad Elements Plot');
%% Plot boundary
% plot(B1(:,1), B1(:,2)), hold on;
% plot(B2(:,1), B2(:,2));
% plot(B3(:,1), B3(:,2));
% plot(B4(:,1), B4(:,2)); 
% legend(["B1","B2","B3","B4"]);
% hold off;
%% Plot solution
% figure
% for i=1:Nels
% edofs=[connectivityData(i,:)];
% fill3 (coordx(edofs),coordy(edofs),u(edofs),u(edofs));hold on
% plot(coordx(edofs),coordy(edofs),'r');hold on
% end
% plot(coordx,coordy,'ro');
% hold off;