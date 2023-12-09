function [Kg, fg] = assembleGlobalMatrixAndForce(coordout, connectivityData)
    Nels = size(connectivityData, 1);
    Nnds = length(coordout(:, 2));
    Kg = zeros(Nnds, Nnds);
    fg = zeros(Nnds, 1);
    switch size(connectivityData, 2)
        case 4
            for i = 1:Nels
                edofs = connectivityData(i, :);
                XN(1:4, 1) = coordout(edofs, 2);
                XN(1:4, 2) = coordout(edofs, 3);
                [Ke, fe] = Elem_Quad4(XN, 0);
                Kg(edofs, edofs) = Kg(edofs, edofs) + Ke;
                fg(edofs, 1) = fg(edofs, 1) + fe;
            end
        case 8
            for i = 1:Nels
                edofs = connectivityData(i, :);
                XN(1:8, 1) = coordout(edofs, 2);
                XN(1:8, 2) = coordout(edofs, 3);
                [Ke, fe] = Elem_Quad8(XN, 0);
                Kg(edofs, edofs) = Kg(edofs, edofs) + Ke;
                fg(edofs, 1) = fg(edofs, 1) + fe;
            end
        otherwise
            error('Unsupported connectivityData size.');
    end
end