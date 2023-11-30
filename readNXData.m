function [coordout, connectivityData] = readNXData(elementsFileName, nodesFileName)
    %% Extract connectivity from NX txt
    fileID0 = fopen(elementsFileName, 'r');
    formatSpec = '%c'; '%d'; %repetido 
    elementLines = splitlines(fscanf(fileID0, formatSpec));
    fclose(fileID0);
    connectivityData = [];
    for i = 1:1:length(elementLines)
        columns = strsplit(elementLines{i});
        cquad4Line = find(contains(columns, 'CQUAD4', 'IgnoreCase', true));
        if ~isempty(cquad4Line)
            numbers = str2double(columns(cquad4Line + 7: cquad4Line + 10));
            connectivityData = [connectivityData; numbers];
        end
    end
    %% Extract nodes from NX txt
    fileID1 = fopen(nodesFileName, 'r');
    formatSpec = '%c'; '%d';
    nos = splitlines(fscanf(fileID1, formatSpec));
    fclose(fileID1);
    aux = [];
    for i= 11:6:length(nos)
        aux = [aux,nos(i)];
    end
    aux = aux';
    colunas = split(aux);
    coord1=[]; % Incialização
    coord2=[];
    coordx=[];
    coordy=[];
    coordout=[];
    for i=1:1:length(colunas)
        coord1=[coord1;colunas(i,5)];
        coordx=[coordx;str2double(coord1(i))];
        coord2=[coord2;colunas(i,6)];
        coordy=[coordy;str2double(coord2(i))];
        coordout=[coordout;i,coordx(i),coordy(i)];
    end
end