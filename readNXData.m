function [coordout, connectivityData] = readNXData(elementsFileName, nodesFileName)
    connectivityData = readElements(elementsFileName);
    coordout = readNodes(nodesFileName);
end
function connectivityData = readElements(elementsFileName)
    % Read connectivity data from elements file
    fileID = fopen(elementsFileName, 'r');
    elementLines = splitlines(fscanf(fileID, '%c')); % Read entire file as a string
    fclose(fileID);
    connectivityData = [];
    
    for i = 1:length(elementLines)
        if contains(elementLines{i}, 'CQUAD4', 'IgnoreCase', true)
            numbers = str2double(strsplit(elementLines{i}));
            elementData = numbers(end-3:end); % Extract the last 4 values for CQUAD4
            connectivityData = [connectivityData; elementData];
        elseif contains(elementLines{i}, 'CQUAD8', 'IgnoreCase', true)
            numbers = str2double(strsplit(elementLines{i}));
            elementData = numbers(end-7:end); % Extract the last 8 values for CQUAD8
            connectivityData = [connectivityData; elementData];
        end
    end
end
function coordout = readNodes(nodesFileName)
    % Read nodes data from nodes file
    fileID = fopen(nodesFileName, 'r');
    nodeLines = splitlines(fscanf(fileID, '%c')); % Read entire file as a string
    fclose(fileID);
    colunas = split(nodeLines(11:6:end));
    coordx_mm = str2double(colunas(:, 5));
    coordy_mm = str2double(colunas(:, 6));
    % Convert from millimeters to meters
    coordx = coordx_mm / 1000;
    coordy = coordy_mm / 1000;
    coordout = [(1:length(coordx))', coordx, coordy];
end