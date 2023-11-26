%PROGRAMA1
close all;
connectivityData;
coordout;
coordx = coordout(:,2);
coordy = coordout(:,3);

k = boundary(coordx,coordy,0.85);
fronteira = [coordx(k),coordy(k)];
plot(fronteira(:,1),fronteira(:,2));

% Identificar partes da fronteira par poder atribuir diferentes condições


%%%%%%%%%%%%% FICOU AQUI



Nels = size(connectivityData, 1);
Nnds = length(coordx);

figure;
hold on;
for i = 1:Nels
    currentConnectivity = connectivityData(i, :);
    x = coordx(currentConnectivity);
    y = coordy(currentConnectivity);
    patch(x, y, 'b', 'FaceAlpha', 0.5);
end
hold off;

axis equal;
grid on;
xlabel('X-axis');
ylabel('Y-axis');
title('Quad Elements Plot');

% inicializacao
Kg=zeros(Nnds,Nnds);
fg=zeros(Nnds,1);

% Carregamento no elemento
for i=1:Nels
    edofs=[connectivityData(i, :)]; %   conectividade deste quad
    % coordenadas do elemento aqui
    XN(1:4,1)=coordx(edofs);
    XN(1:4,2)=coordy(edofs);
    % calculos no elemento
    fL= 4.0;
    [Ke, fe]=Elem_Quad4(XN,fL);
    % assemblagem
    Kg(edofs,edofs)= Kg(edofs,edofs) + Ke;
    fg(edofs,1)= fg(edofs,1) + fe;
end

% Condições de fronteira
boom=1.0e+14;
for i=1:length(k)
    j = k(i); % j é o numero do nó na fronteira
    Kg(j,j) = boom;
    fg(j,j)= boom*0;
end

% Resolução do sistema
u=Kg\fg;

% Boneco
figure
for i=1:Nels
edofs=[connectivityData(i,:)];
fill3 (coordx(edofs),coordy(edofs),u(edofs),u(edofs));hold on
plot(coordx(edofs),coordy(edofs),'r');hold on
end
plot(coordx,coordy,'ro');




