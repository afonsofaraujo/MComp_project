function [Kg, fg] = applyBCs(Kg, fg, essentialBCs)
    boom = 1.0e+14;
    for i = 1:size(essentialBCs, 1)
        Kg(essentialBCs(i,1), essentialBCs(i,1)) = boom;
        fg(essentialBCs(i,1)) = boom*essentialBCs(i,2);
    end
end