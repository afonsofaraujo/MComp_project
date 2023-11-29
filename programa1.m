%PROGRAMA1
close all; % clear; clc
connectivityData;
coordout;
coordx = coordout(:,2);
coordy = coordout(:,3);

k = boundary(coordx,coordy,0.85);
fronteira = [coordx(k),coordy(k)];

% Identificar partes da fronteira para poder atribuir diferentes condições
ymax = max(fronteira(:,2));
xmax = max(fronteira(:,1));
xmin = min(fronteira(:,1));
B1 = [];
B2 = [];
B3 = [];
B4 = [];
s = [];
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
s = unique(s);
for i = 1:size(fronteira,1)
    if ~ismember(i,s)
        B2 = [B2; fronteira(i,:)];
    end
end
B2 = [B2; 3800 600];    % adicionar ponto final
B2 = [0 0; B2];         % adicionar ponto inicial

%%%%%%%%%%%%%

Nels = size(connectivityData, 1);
Nnds = length(coordx);

% inicializacao
Kg=zeros(Nnds,Nnds);
fg=zeros(Nnds,1);
U = 2.5*1000;
boom=1.0e+14;

% Carregamento no elemento
for i=1:Nels
    edofs=[connectivityData(i, :)];
    XN(1:4,1)=coordx(edofs);
    XN(1:4,2)=coordy(edofs);
    fL= 0;
    [Ke, fe]=Elem_Quad4(XN,fL);
    Kg(edofs,edofs)= Kg(edofs,edofs) + Ke;
    fg(edofs,1)= fg(edofs,1) + fe;
end

% B1
for i = 1:size(B1,1)
    j = find(coordx == B1(i, 1) & coordy == B1(i, 2));
    Kg(j, j) = boom;
    fg(j) = boom*coordy(j)*U; % stream function = yU
end
% B2
for i = 1:size(B2,1)
    j = find(coordx == B2(i, 1) & coordy == B2(i, 2)); % Find the corresponding node
    Kg(j, j) = boom;
    fg(j) = 0; % stream function = 0
end
% B4
for i = 1:size(B4,1)
    j = find(coordx == B4(i, 1) & coordy == B4(i, 2));
    Kg(j, j) = boom;
    fg(j) = boom*ymax*U; % stream function = y_maxU
end

% Resolução do sistema
u=Kg\fg;

% % Boneco
% figure
% for i=1:Nels
% edofs=[connectivityData(i,:)];
% fill3 (coordx(edofs),coordy(edofs),u(edofs),u(edofs));hold on
% plot(coordx(edofs),coordy(edofs),'r');hold on
% end
% plot(coordx,coordy,'ro');
% hold off;

%%%%%%%%%%%%%%%%%% pos processamento

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



figure
Res = 0;

% Initialize arrays to store results
arrow_x = [];
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
        
        % Store results for quiver plot
        arrow_x = [arrow_x; xpint(1)];
        arrow_y = [arrow_y; xpint(2)];
        arrow_u = [arrow_u; fluxu(1)];
        arrow_v = [arrow_v; fluxu(2)];
    end
end

% Plot the vector field using quiver
quiver(arrow_x, arrow_y, arrow_u, arrow_v, 1, 'b');  % Set color to blue
hold on;

% Plot other elements
for i = 1:Nels
    edofs = connectivityData(i,:);
    fill(coordx(edofs), coordy(edofs), u(edofs));
    plot(coordx(edofs), coordy(edofs), 'k');
end

hold off;

axis equal;
grid on;
xlabel('X');
ylabel('Y');
title('Vector Field Visualization');


% for i = 1:Nels
%     edofs = connectivityData(i,:);
%     XN(1:4,1) = coordx(edofs);
%     XN(1:4,2) = coordy(edofs);
%     csi = 0;
%     eta = 0;
%     nip = 4;
%     [xp, wp] = Genip2DQ(nip);
% 
%     for ip = 1:nip
%         csi = xp(ip,1);
%         eta = xp(ip,2);
%         [B, psi, Detj] = Shape_N_Der4(XN, csi, eta);
%         uint = psi' * u(edofs);
%         Res = Res + uint * Detj * wp(ip);
%         xpint = XN' * psi;
%         gradu = B' * u(edofs);
%         fluxu = -gradu*0.1;
%         %plot(coordx(edofs), coordy(edofs), 'k');
%         %plot(xpint(1), xpint(2), 'bx');
%         quiver(xpint(1), xpint(2), fluxu(1), fluxu(2), 1/10, 'b');  % Set color to black
%     end
% end
% 










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


%% Plot elementos

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
% plot(fronteira(:,1),fronteira(:,2));

% plot(B1(:,1), B1(:,2)), hold on;
% plot(B2(:,1), B2(:,2));
% plot(B3(:,1), B3(:,2));
% plot(B4(:,1), B4(:,2)); 
% legend(["B1","B2","B3","B4"]);
% hold off;