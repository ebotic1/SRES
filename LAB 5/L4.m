addpath("./funkcije");
defKonstante();

eg = 115e3 / sqrt(3);
f = 50;

bazniNapon1 = 110e3/sqrt(3);
bazniNapon2 = 20e3/sqrt(3);
bazniNapon3 = 0.4e3/sqrt(3);
baznaSnaga = 1e6;

Zb1 = bazniNapon1 * bazniNapon1 / baznaSnaga;
Zb2 = bazniNapon2 * bazniNapon2 / baznaSnaga;
Zb3 = bazniNapon3 * bazniNapon3 / baznaSnaga;

Ib1 = baznaSnaga/bazniNapon1;
Ib2 = baznaSnaga/bazniNapon2;
Ib3 = baznaSnaga/bazniNapon3;

eg=eg/bazniNapon1;

Zg = (0.1 + j*0.5)/Zb1;
Yg = inv(Zg);
Yg = eye(3)*Yg;
I = Yg*(eye(3).*eg);
I = [I(1,1); I(2,2); I(3,3)];

Z0=(0.2095 + j*0.3983)/Zb2;
Y0=(j*2*pi*f*0.0015e-6)*Zb2;
Z1 =(0.0754+j*0.1518)/Zb2;
Y1=(j*2*pi*f*0.002e-6)*Zb2;

linija1 = PiMatricaABC(6.5*Z0, 6.5*Y0, 6.5*Z1, 6.5*Y1);
linija2 = PiMatricaABC(12*Z0, 12*Y0, 12*Z1, 12*Y1);
linija3 = PiMatricaABC(9*Z0, 9*Y0, 9*Z1, 9*Y1);

L1_11 = linija1([1,2,3], 1:3);
L1_12 = linija1([1,2,3], 4:6);
L1_21 = linija1([4,5,6], 1:3);
L1_22 = linija1([4,5,6], 4:6);

L2_11 = linija2([1,2,3], 1:3);
L2_12 = linija2([1,2,3], 4:6);
L2_21 = linija2([4,5,6], 1:3);
L2_22 = linija2([4,5,6], 4:6);

L3_11 = linija3([1,2,3], 1:3);
L3_12 = linija3([1,2,3], 4:6);
L3_21 = linija3([4,5,6], 1:3);
L3_22 = linija3([4,5,6], 4:6);

%transformator 1

Zs0 = izracunajSerijuTransformator(10e6, 110e3, 13.2, 2750) / Zb1;
Ysh0 = izracunajParaleluTransformator(10e6, 110e3, 0.083, 1833) * Zb1;
Zs1 = izracunajSerijuTransformator(10e6, 110e3, 12, 3.3e3) / Zb1;
Ysh1 = izracunajParaleluTransformator(10e6,110e3, 0.2, 2640) * Zb1;
zgi = (50 + j*20)/Zb1;

n = (110/20) * exp(j*5*pi/6);
[Pfin1, Nfin1] = izracunajPozNegTransformatora(Zs1, Ysh1, n);


Zfin1 = [1/(3*zgi+Zs0), 0; 0, 0];
Zfin1 = [1, n; conj(n), n*n] * Zfin1;

T1Y11 = inv(T) * diag([Zfin1(1,1) Pfin1(1,1) Nfin1(1,1)]) * T;
T1Y12 = inv(T) * diag([Zfin1(1,2) Pfin1(1,2) Nfin1(1,2)]) * T;
T1Y21 = inv(T) * diag([Zfin1(2,1) Pfin1(2,1) Nfin1(2,1)]) * T;
T1Y22 = inv(T) * diag([Zfin1(2,2) Pfin1(2,2) Nfin1(2,2)]) * T;

% transformator 2

Zs0 = izracunajSerijuTransformator(1e6, 20e3, 13.2, 0.073e3) / Zb2;
Ysh0 = izracunajParaleluTransformator(1e6, 20e3, 0.083, 0.049e3) * Zb2;
Zs1 = izracunajSerijuTransformator(1e6, 20e3, 12, 0.088e3) / Zb2;
Ysh1 = izracunajParaleluTransformator(1e6,20e3, 0.1, 0.07e3) * Zb2;
zgi = (1+j*0.2)/Zb3;

n = (20/0.4) * exp(j*7*pi/6);
[Pfin2, Nfin2] = izracunajPozNegTransformatora(Zs1, Ysh1, n);


Zfin2 = [0, 0; 0, 1/(3*zgi+Zs0)];
Zfin2 = [1, n; conj(n), n*n] * Zfin2;

T2Y11 = inv(T) * diag([Zfin2(1,1) Pfin2(1,1) Nfin2(1,1)]) * T;
T2Y12 = inv(T) * diag([Zfin2(1,2) Pfin2(1,2) Nfin2(1,2)]) * T;
T2Y21 = inv(T) * diag([Zfin2(2,1) Pfin2(2,1) Nfin2(2,1)]) * T;
T2Y22 = inv(T) * diag([Zfin2(2,2) Pfin2(2,2) Nfin2(2,2)]) * T;

Yc = (j*2*pi*f*200e-9)*Zb2;
Zc = 1/Yc;

Delta = [2*Yc, -Yc, -Yc; -Yc, 2*Yc, -Yc; -Yc, -Yc, 2*Yc];
Zvijezda = [Zc, 0, 0; 0, Zc, 0; 0, 0, Zc];


Ip3 = [188.717373*exp(1i*(-171.65)*pi/180); 
       213.019779*exp(1i*68.99*pi/180); 
       203.895107*exp(1i*(-57.23)*pi/180)];

Ip3 = Ip3./Ib2;

Ip5 = [372.659649*exp(1i*(-23.06)*pi/180);
       425.343276*exp(1i*(-142.50)*pi/180);
       404.906371*exp(1i*90.77*pi/180)];

Ip5 = Ip5./Ib3;


Yrj = [Yg + T1Y11, T1Y12, zeros(3,3), zeros(3,3), zeros(3,3);
       T1Y21, T1Y22 + L1_11+L2_11, L1_12, L2_12, zeros(3,3);
       zeros(3,3), L1_21, L1_22 + Delta + L3_11, L3_12, zeros(3,3);
       zeros(3,3), L2_21, L3_21, L2_22+L3_22+Zvijezda+T2Y11, T2Y12;
       zeros(3,3),zeros(3,3),zeros(3,3), T2Y21, T2Y22];
   
Irj = [I; zeros(3,1); -Ip3; zeros(3,1); -Ip5];

Urj = inv(Yrj) * Irj;












