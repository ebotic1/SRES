function [P, N] = izracunajPozNegTransformatora(Zs1, Ysh1, n)

p = [1/Zs1 + Ysh1/2, -1/Zs1; -1/Zs1, 1/Zs1 + Ysh1/2];

P = [1, 0; 0, conj(n)] * p * [1 , 0; 0, n];
N = conj([1, 0; 0, conj(n)])  * p * conj([1 , 0; 0, n]);
end

