%PROGRAMA3
close all; clear; clc
%[coordout, connectivityData] = readNXData('Elementscasosimples.txt', 'Nodescasosimples.txt');
%[coordout, connectivityData] = readNXData('Elementscomsimetria.txt', 'Nodescomsimetria.txt');
fileName = 'dados-escalar.txt';
[nodeCoordinates, matrixIncidences, materialProperties,... 
distributedLoads, essentialBCs, pointLoads, imposedFlux,...
naturalConvection] = readDadosEscalar(fileName);

coordx = nodeCoordinates(:,2);
coordy = nodeCoordinates(:,3);
connectivityData = matrixIncidences(:, 4:7); % No futuro identificar aqui se é quad ou 8
disp('Data loaded...');

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

% Essential
boom = 1.0e+14;
for i = 1:size(essentialBCs, 1)
    Kg(essentialBCs(i,1), essentialBCs(i,1)) = boom;
    fg(essentialBCs(i,1)) = boom*essentialBCs(i,2);
end
% Neuman
for i = 1:size(imposedFlux,1)
    gama = imposedFlux(i,4);
    fg(imposedFlux(i,2)) = fg(imposedFlux(i,2)) + gama/2;
    fg(imposedFlux(i,3)) = fg(imposedFlux(i,3)) + gama/2;
end
% Robin
for i = 1:size(naturalConvection,1)
    p = naturalConvection(i,4);
    Tout = naturalConvection(i,5);
    Kg(naturalConvection(i,2),naturalConvection(i,2)) = Kg(naturalConvection(i,2),naturalConvection(i,2)) + p/3;
    Kg(naturalConvection(i,2),naturalConvection(i,3)) = Kg(naturalConvection(i,2),naturalConvection(i,3)) + p/6;
    Kg(naturalConvection(i,3),naturalConvection(i,2)) = Kg(naturalConvection(i,3),naturalConvection(i,2)) + p/6;
    Kg(naturalConvection(i,3),naturalConvection(i,3)) = Kg(naturalConvection(i,3),naturalConvection(i,3)) + p/3;
    fg(naturalConvection(i,2)) = fg(naturalConvection(i,2)) + p * Tout / 2;
    fg(naturalConvection(i,3)) = fg(naturalConvection(i,3)) + p * Tout / 2;
end

disp('Solution...');
u=Kg\fg;    % Resolução do sistema

% Plot solution
figure
for i=1:Nels
edofs=[connectivityData(i,:)];
fill3 (coordx(edofs),coordy(edofs),u(edofs),u(edofs));hold on
plot(coordx(edofs),coordy(edofs),'r');hold on
end
plot(coordx,coordy,'ro');
hold off;

%% Plot elements
figure;
hold on;
for i = 1:Nels
    currentConnectivity = connectivityData(i, :);
    x = coordx(currentConnectivity);
    y = coordy(currentConnectivity);
    patch(x, y, 'b', 'FaceAlpha', 0.5);
end
hold off;
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