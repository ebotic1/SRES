function [Ysh] = izracunajParaleluTransformator(Sn, Un, i0, Pfe)

y = (i0*Sn)/(100*Un*Un);
g = Pfe/(Un*Un);
b = sqrt(y*y - g*g);
Ysh = g - j*b;

end

