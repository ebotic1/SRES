function [Z] = izracunajSerijuTransformator(Sn, Un, Uk, Pcu)

z = (Uk*Un*Un)/(100*Sn);
r = (Un*Un*Pcu)/(Sn*Sn);
x = sqrt(z*z - r*r);

Z = r + j*x;

end

