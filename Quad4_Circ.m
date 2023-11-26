clear all
%--------------------------------------------------------------------------
%   Putting things together : Primeiro teste completo dos quads de 4-nos
%   num circulo de raio 1 e representacao grafica
%   Autor: Prof. J. L. M. Fernandes
%--------------------------------------------------------------------------
%
%   definicao "manual" da malha: teste com 3x4=12 quads
quad=[1 2 10 9;2 3 11 10;3 4 12 11;4 5 13 12;5 6 14 13;6 7 15 14;7 8 16 15;
    8 1 9 16;9 10 11 17;11 12 13 17;13 14 15 17;15 16 9 17]
Nquads=size(quad,1)

G=0.5
H=0.6
x1=[1 0.707107 0 -0.707107 -1 -0.707107 0 0.707107 H G 0 -G ...
    -H -G 0 G 0.0]'
y1=[0 0.707107 1 0.707107 0 -0.707107 -1 -0.707107 0 G H G ...
    0 -G -H -G 0.0]'
Nnds=size(x1,1)
%   fim da geracao de malha aqui
%--------------------------------------------------------------------------
%   Tarefa 39, 1a parte : seccao para desenhar a malha de quads-4
%   com visualizacao a 2D da solucao exacta
%--------------------------------------------------------------------------
%   
u=1-x1.^2-y1.^2
figure(1)

Nquads=size(quad,1);
for i=1:Nquads;

edofs=[quad(i,:)]; %   conectividade deste quad
%   Sugestao: por e tirar comentario aqui para se ver a diferenca
%   com fill e sem fill
fill (x1(edofs),y1(edofs),u(edofs));hold on

plot(x1(edofs),y1(edofs),'k');hold on
end
plot(x1,y1,'ro');
%   fim da visualizacao da malha com a solucao exacta
%-----------------------------------------------------
%
Nels=size(quad,1);       % numero de elementos quad-4

Nnds =size(x1,1);        % numero de nos

%-----------------------------------------------------
%   Tarefa 36 - assemblagem dos elementos
%-----------------------------------------------------
%   inicializacao a zeros
Kg=zeros(Nnds,Nnds);
fg=zeros(Nnds,1);
%   
for i=1:Nels   % ciclo para os elementos

edofs=[quad(i,:)]; %   conectividade deste quad
  %     coordenadas do elemento aqui
  XN(1:4,1)=x1(edofs);
  XN(1:4,2)=y1(edofs);
  %     calculos no elemento
fL= 4.0;
[Ke fe]=Elem_Quad4 (XN,fL);

  %     assemblagem
  Kg(edofs,edofs)= Kg(edofs,edofs) + Ke;  % 
  fg(edofs,1)= fg(edofs,1) + fe;          % 
 
end %for i
%Kg
%fg
%-------------------------------------------------------
%   Tarefa 37 - condicoes de fronteira essenciais aqui
%-------------------------------------------------------
boom=1.0e+14;
kount=0;
for i=1:Nnds;
    % calcular o raio
    r =sqrt(x1(i)^2+y1(i)^2);
    if (r > 0.985)  % detecao dos nos da fronteira feita pelo valor do raio
        %   este no esta na fronteira
    Kg(i,i) = boom;
    fg(i,1)= boom*0;
    kount=kount+1;
    end
end
%-----------------------------------------------------
%   Tarefa 38 - seccao para resolver o sistema e
%   controlar o valor maximo
%-----------------------------------------------------
u=Kg\fg;
kount       % numero de nos detectados na fronteira
umx=max(u)
%----------------------------------------------------------
%   Tarefa 39 - seccao para desenhar a solucao aprox. a 3D
%----------------------------------------------------------
figure
%
for i=1:Nels
edofs=[quad(i,:)];
fill3 (x1(edofs),y1(edofs),u(edofs),u(edofs));hold on
plot(x1(edofs),y1(edofs),'r');hold on
end
plot(x1,y1,'ro');
%-----------------------------------------------------------

%-----------------------------------------------------------
%	Tarefa 39 : Visualizacao do erro a 3D usando fill3
%-----------------------------------------------------------
erru=zeros(Nnds,1);	%   assegurar que e 1 vector coluna
      % solucao exacta
      uex=1-x1.^2-y1.^2 ;
      figure
erru = abs(u - uex); % erro com ou sem sinal aqui

for i=1:Nels;
edofs=[quad(i,:)];
fill3 (x1(edofs),y1(edofs),erru(edofs),erru(edofs));hold on
plot(x1(edofs),y1(edofs),'k');hold on
end
plot(x1,y1,'ro');
%-------------- fim da Tarefa 39 -------------------------------------

%---------------------------------------------------------------------
%   Tarefa 40 - calcular (gradiente) fluxo nos pontos de integracao
%---------------------------------------------------------------------

figure
Res=0
for i=1:Nels;
edofs=[quad(i,:)]; %   conectividade deste quad
  XN(1:4,1)=x1(edofs);
  XN(1:4,2)=y1(edofs);
%   O centroide esta na origem
csi=0;
eta=0;
%--------------------------------------------
%   gerar pontos de integracao
nip = 4;    %   experimentar com 9 pontos e vizualizar

[xp wp]=Genip2DQ (nip);

%   percorrer os pontos de integracao
for ip=1:nip;

csi = xp(ip,1);
eta = xp(ip,2);
%--------------------------------------------
%   para cada ponto calcular
%----------------------------------------------------------------
[B psi Detj]=Shape_N_Der4 (XN,csi,eta);
%----------------------------------------------------------------
uint = psi'*u(edofs);

Res = Res + uint*Detj*wp(ip);
xpint = XN'*psi;    %   coordenadas (x,y) do ponto de integracao
gradu = B'*u(edofs);
fluxu = -gradu;
fill (x1(edofs),y1(edofs),u(edofs));hold on
plot(x1(edofs),y1(edofs),'k');hold on
plot(xpint(1),xpint(2),'bx');hold on
quiver(xpint(1),xpint(2),fluxu(1),fluxu(2),1/10);hold on
end

end
Res;
plot(x1,y1,'ro');

errumx=max(erru)
erms= rms(erru)