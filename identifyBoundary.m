function [fronteira, B1, B2, B3, B4] = identifyBoundary(coordout)
    coordx = coordout(:, 2);
    coordy = coordout(:, 3);
    k = boundary(coordx, coordy, 0.85);
    fronteira = [coordx(k), coordy(k)];

    ymax = max(fronteira(:,2));
    xmax = max(fronteira(:,1));
    xmin = min(fronteira(:,1));
    B1 = [];    % esquerda
    B2 = [];    % baixo
    B3 = [];    % direita
    B4 = [];    % cima
    s = [];     % índices já utilizados
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
    s = unique(s);  % retirar índices repetidos (cantos)
    for i = 1:size(fronteira,1)
        if ~ismember(i,s)
            B2 = [B2; fronteira(i,:)]; % B2 contém os pontos da fronteira não utilizados
        end
    end
    B2 = [B2; xmax min(B3(:,2))];   % adicionar ponto final 
    B2 = [0 0; B2];                 % adicionar ponto inicial

    disp('Boundary nodes identified...');
end