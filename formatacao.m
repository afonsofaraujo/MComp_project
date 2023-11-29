%FORMATACAO
format long
clear all
%% Extract connectivity from NX txt
fileID0 = fopen('Elementscomsimetria.txt','r');
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
fileID1 = fopen('Nodescomsimetria.txt','r');
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
%%
% fileID2 = fopen('result.txt','wt');
% writematrix("Escoamento Potencial", 'result.txt','WriteMode', 'append');
% writematrix("Coordenadas dos Nós", 'result.txt','WriteMode', 'append');
% writematrix(length(coordout), 'result.txt','Delimiter', ' ','WriteMode', 'append');
% writematrix(coordout, 'result.txt','Delimiter', ' ','WriteMode', 'append');
% writematrix("Conectividades", 'result.txt','WriteMode', 'append');
% fclose(fileID2);  % Close the file you opened for writing
% 


%% Verificar boundary se trocar o ficheiro
% plot(coordout(:,2),coordout(:,3));
% k=boundary(coordout(:,2),coordout(:,3),0.85);
% fronteira=[coordx(k),coordy(k)];
% plot(fronteira(:,1),fronteira(:,2))

disp('dados carregados')
