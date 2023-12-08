%Função de forma para elementos lineares

function [B, psi, Detj]=Shape_N_Der4 (XN,csi,eta)
    psi=zeros(4,1);
    psi(1) = (1-csi)*(1-eta)/4;
    psi(2) = (1+csi)*(1-eta)/4;
    psi(3) = (1+csi)*(1+eta)/4;
    psi(4) = (1-csi)*(1+eta)/4;
    %derivada de psi em ordem a csi
    Dpsi(1,1) = (eta-1)/4;
    Dpsi(2,1) = (1-eta)/4;
    Dpsi(3,1) = (1+eta)/4;
    Dpsi(4,1) = -(1+eta)/4;
    %derivada de psi em ordem a eta
    Dpsi(1,2) = (csi-1)/4;
    Dpsi(2,2) = -(1+csi)/4;
    Dpsi(3,2) = (1+csi)/4;
    Dpsi(4,2) = (1-csi)/4;
    %jacobiano de csi eta para x y
    jaco = XN'*Dpsi ;
    Detj = det(jaco) ;
    %jacobiano de x y para csi eta
    Invj = inv(jaco) ;
    %matriz com as derivadas das funções de forma em x y
    %(1ª coluna) e em ordem a y (2ª coluna)
    B = Dpsi*Invj ;
end