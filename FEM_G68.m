% FEM_G68
close all;
clear, clc;

%%%%%%%%%%%%%
% inputFileName = 'input_Q4simples.txt';
inputFileName = 'input_Q4base.txt';
% inputFileName = 'input_Q8base.txt';
%%%%%%%%%%%%%

[nodeCoordinates, matrixIncidences, materialProperties,... 
 distributedLoads, essentialBCs, pointLoads, imposedFlux,...
 naturalConvection, elementType, boundaryParameter] = readDadosEscalar(inputFileName);

disp('Load Data...');
coordx = nodeCoordinates(:,2);
coordy = nodeCoordinates(:,3);

if strcmp(elementType, 'QUAD4')
    connectivityData = matrixIncidences(:, 4:7);
elseif strcmp(elementType, 'QUAD8')
    connectivityData = matrixIncidences(:, 4:11);
end

disp('Assembly...');
[Kg, fg] = assemblyGlobalMatrixAndForce(nodeCoordinates, connectivityData);

disp('Boundary Conditions...');
[Kg, fg] = applyBCs(Kg, fg, essentialBCs);

disp('Solution...');
u=Kg\fg;    % Resolução do sistema

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Post-processing...');

[fronteira, B1, B2, B3, B4] = identifyBoundary(nodeCoordinates, boundaryParameter);
[xcentroid, ycentroid] = calculateCentroids(connectivityData, coordx, coordy, elementType);
pressure = calculatePressure(connectivityData, coordx, coordy, u, elementType, materialProperties(1,2));
[vx, vy] = calculateVelocityAtCentroids(connectivityData, coordx, coordy, u, elementType);
[Res, xint, yint, vxint, vyint] = calculateVelocityAtIntegrationPoints(connectivityData, coordx, coordy, u, 4, elementType);

% Output to a text file
randomNumber = randi([10000, 99999]);
outputFileName = ['output_', num2str(randomNumber), '.txt'];
outputFile = fopen(outputFileName, 'w');

% Write data to the file
fprintf(outputFile, 'Título: %s\n', outputFileName);
fprintf(outputFile, 'Tipo de elemento: %s\n', elementType);

% Valores da função corrente nos nós
fprintf(outputFile, 'Valores da função corrente nos nós\n');
fprintf(outputFile, '%d \n', length(u));
for i = 1:length(u)
    fprintf(outputFile, '%f \n', u(i));
end

% Velocidade nos centróides
fprintf(outputFile, 'Velocidades nos centróides\n');
fprintf(outputFile, 'xcentroid, ycentroid, vx, vy\n');
fprintf(outputFile, '%d \n', length(xcentroid));
for i = 1:length(xcentroid)
    fprintf(outputFile, '%f, %f, %f, %f\n', xcentroid(i), ycentroid(i), vx(i), vy(i));
end

% Velocidade nos pontos de integração
fprintf(outputFile, 'Velocidades nos pontos de integração\n');
fprintf(outputFile, 'xint, yint, vxint, vyint\n');
fprintf(outputFile, '%d\n', length(xint));
for i = 1:length(xint)
    fprintf(outputFile, '%f, %f, %f, %f\n', xint(i), yint(i), vxint(i), vyint(i));
end

% Pressão nos centróides
fprintf(outputFile, 'Pressão nos centróides (em bar)\n');
fprintf(outputFile, 'xcentroid, ycentroid, pressure\n');
fprintf(outputFile, '%d\n', length(xcentroid));
for i = 1:length(xcentroid)
    fprintf(outputFile, '%f, %f, %f\n', xcentroid(i), ycentroid(i), pressure(i));
end

fclose(outputFile);
disp(['Results written to ', outputFileName]);
fclose('all');  % Close all open files