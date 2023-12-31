function [B, psi, Detj] = Shape_N_Der8(XN, csi, eta)
    psi = zeros(8,1);
    psi(1) = (csi-1)*(eta+csi+1)*(1-eta)/4;
    psi(2) = (1+csi)*(1-eta)*(csi-eta-1)/4;
    psi(3) = (1+csi)*(1+eta)*(csi+eta-1)/4;
    psi(4) = (csi-1)*(csi-eta+1)*(1+eta)/4;
    psi(5) = (1-csi*csi)*(1-eta)/2;
    psi(6) = (1+csi)*(1-eta*eta)/2;
    psi(7) = (1-csi*csi)*(1+eta)/2;
    psi(8) = (1-csi)*(1-eta*eta)/2;
    Dpsi(1,1) = (2*csi+eta)*(1-eta)/4;
    Dpsi(2,1) = (2*csi-eta)*(1-eta)/4;
    Dpsi(3,1) = (2*csi+eta)*(1+eta)/4;
    Dpsi(4,1) = (2*csi-eta)*(1+eta)/4;
    Dpsi(5,1) = csi*(eta-1);
    Dpsi(6,1) = (1-eta*eta)/2;
    Dpsi(7,1) = -csi*(1+eta);
    Dpsi(8,1) = (eta*eta-1)/2;
    Dpsi(1,2) = (2*eta+csi)*(1-csi)/4;
    Dpsi(2,2) = (2*eta-csi)*(1+csi)/4;
    Dpsi(3,2) = (2*eta+csi)*(1+csi)/4;
    Dpsi(4,2) = (2*eta-csi)*(1-csi)/4;
    Dpsi(5,2) = (csi*csi-1)/2;
    Dpsi(6,2) = -(1+csi)*eta;
    Dpsi(7,2) = (1-csi*csi)/2;
    Dpsi(8,2) = (csi-1)*eta;
    jaco = XN'*Dpsi;
    Detj = det(jaco);
    Invj = inv(jaco);
    B = Dpsi*Invj;
end