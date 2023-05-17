function [Y11, Y12, Y21, Y22] = izracunajPiMatrice(Z0, Y0, Z1, Y1)
Y11 = [(1/Z0)+(Y0/2), 0, 0; 0, (1/Z1)+(Y1/2), 0; 0, 0, (1/Z1) + (Y1/2)];
Y12 = [-1/Z0, 0, 0; 0, -1/Z1, 0; 0, 0, -1/Z1];
Y22 = Y11;
Y21 = Y12;
end

