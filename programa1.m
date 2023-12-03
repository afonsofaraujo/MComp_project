%PROGRAMA1
close all; clear; clc
[coordout, connectivityData] = readNXData('Elementscomsimetria.txt', 'Nodescomsimetria.txt');
disp('Data loaded...');
coordx = coordout(:,2);
coordy = coordout(:,3);
k = boundary(coordx, coordy, 0.85);
fronteira = [coordx(k), coordy(k)];
ymax = max(fronteira(:, 2));
xmax = max(fronteira(:, 1));
xmin = min(fronteira(:, 1));
B4 = fronteira(fronteira(:, 2) == ymax, :);
B1 = fronteira(fronteira(:, 1) == xmin, :);
B3 = fronteira(fronteira(:, 1) == xmax, :);
s = unique([find(fronteira(:, 2) == ymax); find(fronteira(:, 1) == xmin); find(fronteira(:, 1) == xmax)]);
B2 = fronteira(~ismember(1:size(fronteira, 1), s), :);
B2 = [0 0; B2; xmax min(B3(:, 2))];
disp('Assembly...');
Nels = size(connectivityData, 1);
Nnds = length(coordx);
% Kg = sparse(length(coordout),length(coordout)); % inicialização
% fg = sparse(length(coordout),1);
Kg = zeros(Nnds,Nnds); % inicialização
fg = zeros(Nnds,1);
U = 2.5*1000;   % velocidade de entrada em mm/s (porque a malha está em mm)
boom = 1.0e+14;
for i=1:Nels    % Carregamento no elemento
    edofs=[connectivityData(i, :)];
    XN(1:4,1)=coordx(edofs);
    XN(1:4,2)=coordy(edofs);
    [Ke, fe]=Elem_Quad4(XN,0);  % carregamento 0
    Kg(edofs,edofs)= Kg(edofs,edofs) + Ke;
    fg(edofs,1)= fg(edofs,1) + fe;
end
% disp('Boundary conditions...');
% j_B1 = ismember(coordx, B1(:, 1)) & ismember(coordy, B1(:, 2));
% j_B2 = ismember(coordx, B2(:, 1)) & ismember(coordy, B2(:, 2));
% j_B4 = ismember(coordx, B4(:, 1)) & ismember(coordy, B4(:, 2));
% Kg(j_B1, j_B1) = boom;
% fg(j_B1) = boom * coordy(j_B1) * U;
% Kg(j_B2, j_B2) = boom;
% fg(j_B2) = 0;
% Kg(j_B4, j_B4) = boom;
% fg(j_B4) = boom * ymax * U;
%  
disp('Boundary conditions...');
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
figure;
quiver(arrow_x, arrow_y, u2, v2);
hold off;
%%%%%%%%%%%%%


% % Plot contour lines of pressure
% figure;
% subplot(1,2,1);
% quiver(arrow_x, arrow_y, arrow_u, arrow_v); % Plot the vector field using quiver
% hold on;
% contour(XN(:,1), XN(:,2), reshape(pressure, size(coordx)), 'ShowText', 'on');
% hold off;
% axis equal;
% title('u and Pressure Contour');
% 
% % Rotate vectors for the second plot
% subplot(1,2,2);
% quiver(arrow_x, arrow_y, u2, v2);
% hold on;
% contour(XN(:,1), XN(:,2), reshape(pressure, size(coordx)), 'ShowText', 'on');
% hold off;
% axis equal;
% title('u2 and Pressure Contour');












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hold on;
% %Plot other elements
% for i = 1:Nels
%     edofs = connectivityData(i,:);
%     fill(coordx(edofs), coordy(edofs), u(edofs));
%     plot(coordx(edofs), coordy(edofs), 'k');
% end
% 
% hold off;

% figure
% Res=0;
% for i=1:Nels
%     edofs=[connectivityData(i,:)];
%     XN(1:4,1)=coordx(edofs);
%     XN(1:4,2)=coordy(edofs);
%     csi=0;
%     eta=0;
%     nip = 4;
%     [xp, wp]=Genip2DQ (nip);
%     for ip=1:nip
%         csi = xp(ip,1);
%         eta = xp(ip,2);
%         [B, psi, Detj]=Shape_N_Der4 (XN,csi,eta);
%         uint = psi'*u(edofs);
%         Res = Res + uint*Detj*wp(ip);
%         xpint = XN'*psi;
%         gradu = B'*u(edofs);
%         fluxu = -gradu;
%         fill (coordx(edofs),coordy(edofs),u(edofs));hold on
%         plot(coordx(edofs),coordy(edofs),'k');hold on
%         plot(xpint(1),xpint(2),'bx');hold on
%         quiver(xpint(1),xpint(2),fluxu(1),fluxu(2),1/100,'b');hold on
%     end
% end
% Res;
% plot(coordx,coordy,'ro');

% figure
% Res=0;
% for i=1:Nels
%     edofs=[connectivityData(i,:)];
%     XN(1:4,1)=coordx(edofs);
%     XN(1:4,2)=coordy(edofs);
%     csi=0;
%     eta=0;
%     nip = 4;
%     [xp, wp]=Genip2DQ (nip);
%     for ip=1:nip
%         csi = xp(ip,1);
%         eta = xp(ip,2);
%         [B, psi, Detj]=Shape_N_Der4 (XN,csi,eta);
%         uint = psi'*u(edofs);
%         Res = Res + uint*Detj*wp(ip);
%         xpint = XN'*psi;
%         gradu = B'*u(edofs);
%         fluxu = -gradu;
%         fill (coordx(edofs),coordy(edofs),u(edofs));hold on
%         plot(coordx(edofs),coordy(edofs),'k');hold on
%         plot(xpint(1),xpint(2),'bx');hold on
%         quiver(xpint(1),xpint(2),fluxu(1),fluxu(2),1/10);hold on
%     end
% end
% Res;
% plot(coordx,coordy,'ro');

% % Calculate velocity components
% du_dy = gradient(u, coordy);
% dv_dx = -gradient(u, coordx);
% 
% 
% % Interpolate velocity components to the meshgrid points
% U = interp2(coordx, coordy, du_dy, X, Y);
% V = interp2(coordx, coordy, dv_dx, X, Y);
% W = zeros(size(X));  % Z-component is set to zero
% 
% % Plot the velocity vector field
% figure;
% quiver3(X, Y, u, U, V, W, 2, 'LineWidth', 1.5);
% xlabel('X');
% ylabel('Y');
% zlabel('Stream Function (\psi)');
% title('Velocity Vector Field');
% axis equal;
% grid on;
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

% Plot boundary
% plot(fronteira(:,1),fronteira(:,2));
% 
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