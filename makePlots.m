% makePlots
close all;
clear, clc;

%%%%%%%%%%%%%
% inputFileName = 'input_Q4simples.txt';
% outputfileName = 'output_Q4simples.txt';
inputFileName = 'input_Q4base.txt';
outputfileName = 'output_Q4base.txt';
% inputFileName = 'input_Q8base.txt';
% outputfileName = 'output_Q8base.txt';
%%%%%%%%%%%%%

[nodeCoordinates, matrixIncidences, materialProperties,... 
 distributedLoads, essentialBCs, pointLoads, imposedFlux,...
 naturalConvection, elementType, boundaryParameter] = readDadosEscalar(inputFileName);
[u, xcentroid, ycentroid, vx, vy, xint, yint, vxint, vyint, pressure] = readOutput(outputfileName);
[fronteira, B1, B2, B3, B4] = identifyBoundary(nodeCoordinates, boundaryParameter);
coordx = nodeCoordinates(:,2);
coordy = nodeCoordinates(:,3);

half = 0; % For half piece half = 1, whole half = 0,

if half

    % Interpolaçao
    [fineX, fineY, fineU, fineP] = interpolateAndMask(coordx, coordy, u, pressure, xcentroid, ycentroid, fronteira, 2);
    
    % Linhas de corrente
    hFig1 = figure;
    contourf(fineX, fineY, fineU, 20), hold on;
    plot(fronteira(:,1), fronteira(:,2), 'Color', 'k', 'LineWidth', 1);
    title('Função corrente'), xlabel('X-axis'), ylabel('Y-axis'), axis equal, axis off, hold off;
    set(gcf, 'Color', 'w');      % Set background color to white
    set(hFig1, 'Name', 'Função Corrente');
    
    % Pressao
    hFig2 = figure;
    contour(fineX, fineY, fineP, 10,'LineColor','k',"ShowText",true,"LabelFormat","%0.2f",'LabelSpacing', 400), hold on;
    plot(fronteira(:,1), fronteira(:,2), 'Color', 'k', 'LineWidth', 1);
    title('Pressão (bar)'), xlabel('X-axis'), ylabel('Y-axis'), axis equal, axis off, hold off;
    set(gcf, 'Color', 'w');      % Set background color to white
    set(hFig2, 'Name', 'Pressão');
    
    % Velocidade
    hFig3 = figure;
    subplot(2, 1, 1);
    quiver(xint, yint, vxint, vyint);
    title('Pontos de integração'), xlabel('X-axis'), ylabel('Y-axis'), axis equal, axis off;
    subplot(2, 1, 2);
    quiver(xcentroid, ycentroid, vx, vy);
    title('Centróides'), xlabel('X-axis'), ylabel('Y-axis'), axis equal, axis off;
    set(gcf, 'Color', 'w');      % Set background color to white
    set(hFig3, 'Name', 'Velocidade');
else

    % Espelhar pontos e resultados
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
    pxcentroid = [xcentroid; xcentroid];
    pycentroid = [ycentroid; 2*max(coordy)-ycentroid];
    pressuret = [pressure; pressure];
    xintt = [xint; xint];
    yintt = [yint; 2*max(coordy)-yint];
    vxintt = [vxint; vxint];
    vyintt = [vyint; -vyint];
    vxt = [vx; vx];
    vyt = [vy; -vy];

    % Interpolaçao
    [fineXt, fineYt, fineUt, finePt] = interpolateAndMask(px, py, ut, pressuret, pxcentroid, pycentroid, fronteirat, 2);
    
    % Linhas de corrente
    hFig1 = figure; 
    contourf(fineXt, fineYt, fineUt, 20), hold on;
    plot([min(px), max(px)], [0.9, 0.9], 'k-', 'LineWidth', 0.5);
    plot(fronteirat(:,1), fronteirat(:,2), 'Color', 'k', 'LineWidth', 1);
    title('Função Corrente'), xlabel('X-axis'), ylabel('Y-axis'), axis equal, axis off, hold off;
    set(gcf, 'Color', 'w');      % Set background color to white
    set(hFig1, 'Name', 'Função Corrente');
    
    % Pressao
    hFig2 = figure;
    contour(fineXt, fineYt, finePt, 10,'LineColor','k',"ShowText",true,"LabelFormat","%0.2f",'LabelSpacing', 400), hold on;
    plot(fronteirat(:,1), fronteirat(:,2), 'Color', 'k', 'LineWidth', 1);
    title('Pressão (bar)'), xlabel('X-axis'), ylabel('Y-axis'), axis equal, axis off, hold off;
    set(gcf, 'Color', 'w');      % Set background color to white
    set(hFig2, 'Name', 'Pressão');

    % Velocidade
    hFig3 = figure;
    subplot(2, 1, 1);
    quiver(xintt, yintt, vxintt, vyintt), hold on;
    B6 = [B2(:,1), 2*max(coordy)-B2(:,2)];
    plot(B2(:,1), B2(:,2),'k', 'LineWidth', 0.5);
    plot(B6(:,1), B6(:,2),'k', 'LineWidth', 0.5);
    %plot(fronteirat(:,1), fronteirat(:,2), 'Color', 'k', 'LineWidth', 0.5);
    title('Pontos de integração'), xlabel('X-axis'), ylabel('Y-axis'), axis equal, axis off, hold off;
    subplot(2, 1, 2);
    quiver(pxcentroid, pycentroid, vxt, vyt), hold on;
    plot(B2(:,1), B2(:,2),'k', 'LineWidth', 0.5);
    plot(B6(:,1), B6(:,2),'k', 'LineWidth', 0.5);
    %plot(fronteirat(:,1), fronteirat(:,2), 'Color', 'k', 'LineWidth', 0.5);
    title('Centróides'), xlabel('X-axis'), ylabel('Y-axis'), axis equal, axis off, hold off;
    set(gcf, 'Color', 'w');      % Set background color to white
    set(hFig3, 'Name', 'Velocidade');
end