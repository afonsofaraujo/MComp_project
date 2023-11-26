clear all
close all
%--------------------------------------------------------------------------
%   Putting things together : Geracao automatica de quads de 4-nos num
%   circulo de raio 1, obtencao da solucao pelo MEF, calculo do erro e
%   representacao grafica.
%   Geracao Automatica por blocos: define 5 blocos fixos de 8 nos para
%   preencher e refinar o circulo de raio 1 com quadrilateros Quad-4.
%   Autor: Prof. J. L. M. Fernandes
%--------------------------------------------------------------------------
M = 8   % 8 pontos no circulo exterior
Nnds = 2*M + M/2  % numero de nos gerados para definir os blocos = 20
r = 1;
alfa=2*pi/M;
%   assegurar que os vectores inicializados a zero são vect. coluna 
x=zeros(Nnds,1);
y=zeros(Nnds,1);
L = 0;
theta = 0;     %   fiada da fronteira
for i=1:M      %   ciclo para os angulos
    L = L + 1;
    x(L) = r*cos(theta);
    y(L) = r*sin(theta);
    theta = theta + alfa;
end
H=0.6   %   parametros do bloco central
G=0.5   %   parametros do bloco central
x(13)=H
x(14)=G
x(15)=0
x(16)=-G
x(17)=-H
x(18)=-G
x(19)=0
x(20)=G
y(13)=0
y(14)=G
y(15)=H
y(16)=G
y(17)=0
y(18)=-G
y(19)=-H
y(20)=-G
%   nos no meio dos lados
x(9)=(x(2)+x(14))/2
y(9)=(y(2)+y(14))/2
x(10)=(x(4)+x(16))/2
y(10)=(y(4)+y(16))/2
x(11)=(x(6)+x(18))/2
y(11)=(y(6)+y(18))/2
x(12)=(x(8)+x(20))/2
y(12)=(y(8)+y(20))/2
%   criar 5 blocos de Quads-8 para preencher o circulo
quad8=[8 2 14 20 1 9 13 12; 2 4 16 14 3 10 15 9; 4 6 18 16 5 11 17 10;
    6 8 20 18 7 12 19 11; 18 20 14 16 19 13 15 17];

%----------- visualizacao rudimentar dos blocos de 8 aqui -----------
%       nao mostra a curvatura existente nos lados curvos

u=1-x.^2-y.^2;

figure(1);
Nels=size(quad8,1);
for i=1:Nels;
    no1=quad8(i,1);
    no2=quad8(i,2);
    no3=quad8(i,3);
    no4=quad8(i,4);
    no5=quad8(i,5);
    no6=quad8(i,6);
    no7=quad8(i,7);
    no8=quad8(i,8);
edofs=[no1 no5 no2 no6 no3 no7 no4 no8]; % desenhar por esta ordem

fill (x(edofs),y(edofs),u(edofs));hold on
plot(x(edofs),y(edofs),'b');hold on
end
plot(x,y,'ro');
%----------- fim da visualizacao dos blocos de 8 aqui -------------

%----------- inicio da geracao da malha de quads-4 aqui -----------
for i=1:Nels;
    no1=quad8(i,1);
    no2=quad8(i,2);
    no3=quad8(i,3);
    no4=quad8(i,4);
    no5=quad8(i,5);
    no6=quad8(i,6);
    no7=quad8(i,7);
    no8=quad8(i,8);
%   ordem natural : de 1 a 8    
edofs=[no1 no2 no3 no4 no5 no6 no7 no8]; % ordem natural
%----------- especificar numero de divisoes aqui --------------------------

%---------
Nx=6;
Ny=2;
%---------
Nx=3;
Ny=1;

%---------
Nx=6;
Ny=2;
%---------
Nx=12;
Ny=4;


Nx=24;
Ny=8;
%--------------------------------------------------------------------------
if (i ==5)
    Ny=Nx;      %   respeitar a simetria no bloco central
end
%--------------------------------------------------------------------------
if (i > 1)
    [quad2 x2 y2]=Block8_quad4 (x(edofs),y(edofs),Nx,Ny);
    [quad x1 y1]=Block_Merge (x1,y1,quad,x2,y2,quad2);
else
[quad x1 y1]=Block8_quad4 (x(edofs),y(edofs),Nx,Ny);
end
end

Nnds =size(x1,1);        % numero de nos dos quads gerados
%----------- corrigir a posicao dos nos da fronteira aqui -----------------
%            para assegurar ter o raio exactamente r=1
kount=0
for i=1:Nnds;
    xx=x1(i);
    yy=y1(i);
    % calcular o raio
    r =sqrt(xx^2+yy^2);
    if (r > 0.98)  % seleccao dos nos da fronteira feita pelo valor do raio
        %   este no esta na fronteira
        %   corrigir as coordenadas para ter o raio exactamente r=1
        kount=kount+1;
        x1(i)=xx/r;
        y1(i)=yy/r;
    end
end

%----------- fim da geracao da malha de quads-4 aqui -----------
%-----------------------------------------------------
%   Tarefa 39 - seccao para desenhar a malha de quads-4
%   com a solucao exacta a 2D
%-----------------------------------------------------
%   
u=1-x1.^2-y1.^2
figure
Nquads=size(quad,1);
for i=1:Nquads;

edofs=[quad(i,:)]; %   conectividade deste quad
fill (x1(edofs),y1(edofs),u(edofs));hold on
plot(x1(edofs),y1(edofs),'r');hold on
end
plot(x1,y1,'ro');
%   fim da visualizacao da malha com a solucao exacta
%-----------------------------------------------------
%
Nels=size(quad,1);       % numero de elementos quads

Nnds =size(x1,1);        % numero de nos

%-----------------------------------------------------
%   Tarefa 36 - assemblagem dos elementos
%-----------------------------------------------------
%   inicializacao a zeros
Kg=zeros(Nnds,Nnds);
fg=zeros(Nnds,1);
%   
for i=1:Nels;

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
%-----------------------------------------------------
%Kg
%fg
%-----------------------------------------------------
%   Tarefa 37 - condicoes de fronteira essenciais
%-----------------------------------------------------
%   
boom=1.0e+12;
kount=0;
for i=1:Nnds;
    % calcular o raio
    r =sqrt(x1(i)^2+y1(i)^2);
    if (r > 0.985) % seleccao dos nos da fronteira feita pelo valor do raio
        %   este no esta na fronteira
    Kg(i,i) = boom;
    fg(i,1)= boom*0;
    kount=kount+1;
    end
end
%-----------------------------------------------------
%   Tarefa 38 - seccao para resolver o sistema e
%   controlar o valor maximo
%----------------------------------------------------
u=Kg\fg;
kount;
umx=max(u)
%-----------------------------------------------------
%   Tarefa 39 - seccao para desenhar a solucao em 3D
%-----------------------------------------------------
figure
%
for i=1:Nels;
edofs=[quad(i,:)];
fill3 (x1(edofs),y1(edofs),u(edofs),u(edofs));hold on
plot(x1(edofs),y1(edofs),'r');hold on
end
plot(x1,y1,'ro');
%--------------------------------------------------------------------------
%
%-------------------------------------------------------
%   Tarefa 39 : Visualizacao do erro em 3D usando fill3  
%-------------------------------------------------------
%	
erru=zeros(Nnds,1);	%   assegurar que e 1 vector coluna
      % solucao exacta
      uex=1-x1.^2-y1.^2 ;
erru = abs(u - uex); % erro com ou sem sinal aqui
ermx= max(erru)
%  normer = norm(erru)/sqrt(Nnds) % esta e igual a rms
erms= rms(erru)
figure      

for i=1:Nels;
edofs=[quad(i,:)];
fill3 (x1(edofs),y1(edofs),erru(edofs),erru(edofs));hold on
plot(x1(edofs),y1(edofs),'k');hold on
end
plot(x1,y1,'ro');
%-------------- fim da Tarefa 39 ---------------------------

%-----------------------------------------------------------
%   Tarefa 40 - calcular (gradiente) fluxo nos centroides
%-----------------------------------------------------------
figure
for i=1:Nels;
edofs=[quad(i,:)]; %   conectividade deste quad
  XN(1:4,1)=x1(edofs);
  XN(1:4,2)=y1(edofs);
%   O centroide esta na origem
csi=0;
eta=0;
%   para cada centroide, calcular
%----------------------------------------------------------------
[B psi Detj]=Shape_N_Der4 (XN,csi,eta);
%----------------------------------------------------------------
uint = psi'*u(edofs);
xpint = XN'*psi;    %   posicao (x,y) do centroide
gradu = B'*u(edofs);
fluxu = -gradu;
%fill (x1(edofs),y1(edofs),u(edofs));hold on
plot(x1(edofs),y1(edofs),'k');hold on
plot(xpint(1),xpint(2),'bx');hold on
quiver(xpint(1),xpint(2),fluxu(1),fluxu(2),1/10);hold on
end
%plot(x1,y1,'ro');


function [quad x y]=Block8_quad4 (xg,yg,Nx,Ny)
%   preenche um bloco de 8 nos com elementos de 4 nos (quad-4)
%   Autor: Prof. J. L. M. Fernandes
%--------------------------------------
Nels=Nx*Ny;
Nnds=(Nx+1)*(Ny+1);
x=zeros(Nnds,1);
y=zeros(Nnds,1);
%   gera lista de pontos
stepx=2/Nx;
stepy=2/Ny;
L=0;
eta = -1;
for j=1: Ny+1;
csi=-1;
for i=1: Nx+1;
    
psi(1) = -(1-csi)*(1-eta)*(csi+eta+1)/4;
psi(2) = (1+csi)*(1-eta)*(csi-eta-1)/4;
psi(3) = (1+csi)*(1+eta)*(csi+eta-1)/4;
psi(4) = (1-csi)*(1+eta)*(eta-csi-1)/4;
psi(5) = (1-csi*csi)*(1-eta)/2;
psi(6) = (1+csi)*(1-eta*eta)/2;
psi(7) = (1-csi*csi)*(1+eta)/2;
psi(8) = (1-csi)*(1-eta*eta)/2;
L=L+1;
x(L) = psi*xg;
y(L) = psi*yg;
csi= csi + stepx;
end
eta= eta + stepy;
end
L=0;
%   gera lista de quads
for j=1: Ny;
    for i=1: Nx;
        L1= (j-1)*(Nx+1) + i;
        L4= j*(Nx+1) + i;
        L=L+1;
        quad(L,1)= L1;
        quad(L,2)= L1+1;
        quad(L,3)= L4 + 1;
        quad(L,4)= L4;
    end
end
end

function [quad x y]=Block_Merge (x1,y1,quad1,x2,y2,quad2)
%--------------------------------------
%   junta dois blocos de elementos do mesmo tipo, gerados separadamente,
%   numa malha unica. Detecta nos comuns e retem o de menor numero
%   Autor: Prof. J. L. M. Fernandes
%--------------------------------------
%   
Nnds1=size(x1,1);
Nnds2=size(x2,1);
x(1:Nnds1,1)=x1(1:Nnds1);
y(1:Nnds1,1)=y1(1:Nnds1);
Nnds=Nnds1;
for j=1:Nnds2;
for i=1:Nnds1;
    % cacula distancia entre dois nos
    dist = (x2(j)-x1(i))^2 +(y2(j)-y1(i))^2;
    elim=0;
    if (dist < 1.0e-8)
    %   detectados nos comuns
    New(j)=i;
    elim=1;
    break;
    end
end
if (elim == 0)
    Nnds = Nnds+1;
    x(Nnds)=x2(j);
    y(Nnds)=y2(j);
    New(j)=Nnds;
end
end
  %     update quad2 list with new numbers
Nels2=size(quad2,1);
for j=1:Nels2;
no1=quad2(j,:);
quad2(j,:)=New(no1);
end
%   forma lista unica de quads
quad=vertcat(quad1,quad2);
end
