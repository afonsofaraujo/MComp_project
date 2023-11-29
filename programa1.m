%PROGRAMA1
close all;
connectivityData;
coordout;
coordx = coordout(:,2);
coordy = coordout(:,3);

k = boundary(coordx,coordy,0.85);
fronteira = [coordx(k),coordy(k)];
%plot(fronteira(:,1),fronteira(:,2));

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

plot(B1(:,1), B1(:,2)), hold on;
plot(B2(:,1), B2(:,2));
plot(B3(:,1), B3(:,2));
plot(B4(:,1), B4(:,2)); 
legend(["B1","B2","B3","B4"]);
hold off;

%%%%%%%%%%%%%

Nels = size(connectivityData, 1);
Nnds = length(coordx);

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

% inicializacao
Kg=zeros(Nnds,Nnds);
fg=zeros(Nnds,1);
U = 2.5/1000;
boom = 10^12;



% Carregamento no elemento
for i=1:Nels
    edofs=[connectivityData(i, :)]; %   conectividade deste quad
    % coordenadas do elemento aqui
    XN(1:4,1)=coordx(edofs);
    XN(1:4,2)=coordy(edofs);
    % calculos no elemento
    fL= 0;
    [Ke, fe]=Elem_Quad4(XN,fL);
    % assemblagem
    Kg(edofs,edofs)= Kg(edofs,edofs) + Ke;
    fg(edofs,1)= fg(edofs,1) + fe;
end
% B2 stream function is zero
for i = 1:size(B2, 1)
    j = find(coordx == B2(i, 1) & coordy == B2(i, 2)); % Find the corresponding node
    Kg(j, j) = 1*boom;
    fg(j) = 0; % stream function = 0
end
% B1
for i = 1:size(B1)
    j = find(coordx == B1(i, 1) & coordy == B1(i, 2));
    Kg(j, j) = boom*1;
    fg(j) = boom*coordy(j)*U; % stream function = yU
end

% B4  - é igual ao B2 ?
for i = 1:size(B1)
    j = find(coordx == B1(i, 1) & coordy == B1(i, 2));
    Kg(j, j) = boom*1;
    fg(j) = boom*coordy(j)*U; % stream function = yU
end


% B3 derivada da função corrente é 0







% % Carregamento no elemento
% for i=1:Nels
%     edofs=[connectivityData(i, :)]; %   conectividade deste quad
%     % coordenadas do elemento aqui
%     XN(1:4,1)=coordx(edofs);
%     XN(1:4,2)=coordy(edofs);
%     % calculos no elemento
%     fL= 4.0;
%     [Ke, fe]=Elem_Quad4(XN,fL);
%     % assemblagem
%     Kg(edofs,edofs)= Kg(edofs,edofs) + Ke;
%     fg(edofs,1)= fg(edofs,1) + fe;
% end
% 
% % Condições de fronteira
% boom=1.0e+14;
% for i=1:length(k)
%     j = k(i); % j é o numero do nó na fronteira
%     Kg(j,j) = boom;
%     fg(j,j)= boom*0;
% end

% Resolução do sistema
u=Kg\fg;

%% Boneco
figure
for i=1:Nels
edofs=[connectivityData(i,:)];
fill3 (coordx(edofs),coordy(edofs),u(edofs),u(edofs));hold on
plot(coordx(edofs),coordy(edofs),'r');hold on
end
plot(coordx,coordy,'ro');




