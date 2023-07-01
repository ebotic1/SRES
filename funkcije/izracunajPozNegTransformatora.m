function [P, N] = izracunajPozNegTransformatora(Zs1, Ysh1, n)

p = [1/Zs1 + Ysh1/2, -1/Zs1; -1/Zs1, 1/Zs1 + Ysh1/2];

P = [1, n; conj(n), n*n] * p;
N = [1, conj(n); n, n*n] * p;
end

