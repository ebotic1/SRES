function [mat] = PiMatricaABC(Z0, Y0, Z1, Y1)

global T;
[rez1 rez2 rez3 rez4] = izracunajPiMatrice(Z0, Y0, Z1, Y1);
mat = {rez1, rez2, rez3, rez4};
for i=1:length(mat)
    mat{i} = inv(T) * mat{i} * T;
end
mat = [mat{1} mat{2}; mat{3} mat{4}];


end

