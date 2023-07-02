addpath("./../funkcije");
defKonstante();

eg = 115e3/sqrt(3);
f = 50;

bazniNapon1 = 110e3/sqrt(3);
bazniNapon2 = 20e3/sqrt(3);
bazniNapon3 = 0.4e3/sqrt(3);
baznaSnaga = 1e7;

Zb1 = bazniNapon1 * bazniNapon1 / baznaSnaga;
Zb2 = bazniNapon2 * bazniNapon2 / baznaSnaga;
Zb3 = bazniNapon3 * bazniNapon3 / baznaSnaga;

Ib1 = baznaSnaga/bazniNapon1;
Ib2 = baznaSnaga/bazniNapon2;
Ib3 = baznaSnaga/bazniNapon3;

eg=eg/bazniNapon1;
eg = [eg; eg*a*a; eg*a];
Zg = (0.1 + j*0.5)/Zb1;
Yg = inv(eye(3).*Zg);
Ig = Yg*eg;

Z0 =(0.2095 + 0.3983*j)/Zb2;
Y0 =(j*2*pi*f*0.0015e-6)*Zb2;
Z1 =( 0.0754 + 0.1518*j)/Zb2;
Y1 =(j*2*pi*f*0.002e-6)*Zb2;

linija1 = PiMatricaABC(6.5*Z0, 6.5*Y0, 6.5*Z1, 6.5*Y1);
linija2 = PiMatricaABC(12*Z0, 12*Y0, 12*Z1, 12*Y1);
linija3 = PiMatricaABC(9*Z0, 9*Y0, 9*Z1, 9*Y1);

[L1_11, L1_12, L1_21, L1_22] = spakujMatricu(linija1);
[L2_11, L2_12, L2_21, L2_22] = spakujMatricu(linija2);
[L3_11, L3_12, L3_21, L3_22] = spakujMatricu(linija3);


%transformator 1

snaga = 10e6;

Zs0 = izracunajSerijuTransformator(snaga, 110e3, 13.2, 2750) / Zb1;
Ysh0 = izracunajParaleluTransformator(snaga, 110e3, 0.083, 1833) * Zb1;
Zs1 =  izracunajSerijuTransformator(snaga, 110e3, 12, 3300) / Zb1;
Ysh1 = izracunajParaleluTransformator(snaga, 110e3, 0.1, 2640) * Zb1;
zgi = (50 + j*20)/Zb1;

n = 1 * exp(j*5*pi/6);
[Pfin1, Nfin1] = izracunajPozNegTransformatora(Zs1, Ysh1, n);


Zfin1 = [1/(3*zgi+Zs0), 0; 0, 0];
Zfin1 = [1, 0; 0, conj(n)] * Zfin1 * [1 , 0; 0, n];

T1Y11 = inv(T) * diag([Zfin1(1,1) Pfin1(1,1) Nfin1(1,1)]) * T;
T1Y12 = inv(T) * diag([Zfin1(1,2) Pfin1(1,2) Nfin1(1,2)]) * T;
T1Y21 = inv(T) * diag([Zfin1(2,1) Pfin1(2,1) Nfin1(2,1)]) * T;
T1Y22 = inv(T) * diag([Zfin1(2,2) Pfin1(2,2) Nfin1(2,2)]) * T;



%transformator 2

snaga = 1e6;

Zs0 = izracunajSerijuTransformator(snaga, 20e3, 13.2, 73) / Zb2;
Ysh0 = izracunajParaleluTransformator(snaga, 20e3, 0.083, 49) * Zb2;
Zs1 =  izracunajSerijuTransformator(snaga, 20e3, 12, 88) / Zb2;
Ysh1 = izracunajParaleluTransformator(snaga, 20e3, 0.1, 70) * Zb2;
zgi = 1+j*0.2;

n = (1) * exp(j*7*pi/6);
[Pfin2, Nfin2] = izracunajPozNegTransformatora(Zs1, Ysh1, n);


Zfin2 = [0, 0; 0, 1/(3*zgi+Zs0)];
Zfin2 = [1, 0; 0, conj(n)] * Zfin2 * [1 , 0; 0, n];

T2Y11 = inv(T) * diag([Zfin2(1,1) Pfin2(1,1) Nfin2(1,1)]) * T;
T2Y12 = inv(T) * diag([Zfin2(1,2) Pfin2(1,2) Nfin2(1,2)]) * T;
T2Y21 = inv(T) * diag([Zfin2(2,1) Pfin2(2,1) Nfin2(2,1)]) * T;
T2Y22 = inv(T) * diag([Zfin2(2,2) Pfin2(2,2) Nfin2(2,2)]) * T;


Yc = (j*2*pi*f*200e-9)*Zb2;
Delta = [2*Yc, -Yc, -Yc; -Yc, 2*Yc, -Yc; -Yc, -Yc, 2*Yc];
Zvijezda = diag([Yc Yc Yc]);


Ip3 = [188.717373*exp(1i*(-171.65)*pi/180); 
       213.019779*exp(1i*68.99*pi/180); 
       203.895107*exp(1i*(-57.23)*pi/180)];

Ip3 = Ip3./Ib2;

Ip5 = [372.659649*exp(1i*(-23.06)*pi/180);
       425.343276*exp(1i*(-142.50)*pi/180);
       404.906371*exp(1i*90.77*pi/180)];
   
Ip5 = Ip5./Ib3;

Yrj = [T1Y11 + Yg, T1Y12, zeros(3,3), zeros(3,3), zeros(3,3);
       T1Y21, T1Y22 + L1_11 + L2_11, L1_12, L2_12, zeros(3,3);
       zeros(3,3), L1_21, L1_22 + L3_11 + Delta, L3_12, zeros(3,3);
       zeros(3,3), L2_21, L3_21, L2_22 + L3_22 + Zvijezda + T2Y11, T2Y12;
       zeros(3,3), zeros(3,3), zeros(3,3), T2Y21, T2Y22];
   
Irj = [Ig; zeros(3,1); -Ip3; zeros(3,1); -Ip5];

Urj = inv(Yrj) * Irj;









