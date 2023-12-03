function [nodeCoordinates, matrixIncidences, materialProperties,... 
        distributedLoads, essentialBCs, pointLoads, imposedFlux,...
        naturalConvection] = readDadosEscalar(fileName)
    %fileName = 'dados-escalar.txt';
    % Extract node coordinates from txt
    fid = fopen(fileName, 'r');
    % Skip the first three lines
    for i = 1:2
        fgetl(fid);
    end
    % Coordenadas dos nos
    Nnds = fscanf(fid, '%d', 1);
    fgetl(fid);
    nodeCoordinates = fscanf(fid, '%d %f %f', [3, Nnds])';
    fgetl(fid);
    fgetl(fid);
    % Matrix de incidencias/conectividades
    Nels = fscanf(fid, '%d', 1);
    fgetl(fid);
    matrixIncidences = fscanf(fid, '%d', [7, Nels])';
    fgetl(fid);
    fgetl(fid);
    % Propriedades/Material
    NTM = fscanf(fid, '%d', 1);
    fgetl(fid);
    if NTM~=0
        materialProperties = fscanf(fid, '%d %f', [2, NTM])';
    else
        materialProperties = [];
    end
    fgetl(fid);
    fgetl(fid);
    % Fontes/carregamentos distribuídos
    NTECD = fscanf(fid, '%d', 1);
    if NTECD~=0
        distributedLoads = fscanf(fid, '%d %f', [2, NTECD])';
    else
        distributedLoads = [];
    end
    fgetl(fid);
    fgetl(fid);
    % Condição fronteira essencial 
    NTGLI = fscanf(fid, '%d', 1);
    fgetl(fid);
    if NTGLI~=0
        essentialBCs = fscanf(fid, '%d %f', [2, NTGLI])';
    else
        essentialBCs = [];
    end
    fgetl(fid);
    fgetl(fid);
    % Fontes/cargas pontuais impostas  
    NTCPI = fscanf(fid, '%d', 1);
    fgetl(fid);
    if NTCPI~=0
        pointLoads = fscanf(fid, '%d %f', [2, NTCPI])';
    else
        pointLoads = [];
    end
    fgetl(fid);
    fgetl(fid);
    % Fluxo imposto na fronteira
    NTEFI = fscanf(fid, '%d', 1);
    fgetl(fid);
    if NTEFI~=0
        imposedFlux = fscanf(fid, '%d %f', [4, NTEFI])';
    else
        imposedFlux = [];
    end
    fgetl(fid);
    fgetl(fid);  
    % CF convecção natural
    NTECN = fscanf(fid, '%d', 1);
    fgetl(fid);
    if NTECN~=0
        naturalConvection = fscanf(fid, '%d %f', [5, NTECN])';
    else
        naturalConvection = [];
    end
    fclose(fid);
end
