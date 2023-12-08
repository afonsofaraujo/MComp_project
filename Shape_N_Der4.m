%Função de forma para elementos lineares

function [B, psi, Detj]=Shape_N_Der4 (XN,csi,eta)
    psi=zeros(4,1);
    psi(1) = (1-csi)*(1-eta)/4;
    psi(2) = (1+csi)*(1-eta)/4;
    psi(3) = (1+csi)*(1+eta)/4;
    psi(4) = (1-csi)*(1+eta)/4;
    Dpsi(1,1) = (eta-1)/4;
    Dpsi(2,1) = (1-eta)/4;
    Dpsi(3,1) = (1+eta)/4;
    Dpsi(4,1) = -(1+eta)/4;
    Dpsi(1,2) = (csi-1)/4;
    Dpsi(2,2) = -(1+csi)/4;
    Dpsi(3,2) = (1+csi)/4;
    Dpsi(4,2) = (1-csi)/4;
    jaco = XN'*Dpsi ;
    Detj = det(jaco) ;
    Invj = inv(jaco) ;
    B = Dpsi*Invj ;
end