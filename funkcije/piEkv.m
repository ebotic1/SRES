function [Y] = piEkv(Z0, Y0)

Y = [1/Z0 + Y0/2, -1/Z0;
    -1/Z0, 1/Z0 + Y0/2];


end

