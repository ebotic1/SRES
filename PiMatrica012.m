function [matrica] = PiMatrica012(Z0, Y0, Z1, Y1)
[rez1 rez2 rez3 rez4] = izracunajPiMatrice(Z0, Y0, Z1, Y1);
matrica = [rez1, rez2; rez3, rez4];
end

