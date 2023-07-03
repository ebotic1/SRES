addpath("./../funkcije");
defKonstante();


slack = 1.05;


napon = [1.05; 1.05; 1.05; 1.05];

p = [slack * conj(slack); -(0.4 + j*0.2); -(0.4 + j*0.2); 0.9];

I = izracunajStruje(p, napon);

Z0 = 0.003 + j*0.006;

L1 = piEkv(Z0, 0);
L2 = L1;
L3 = L1;


Y = [1, 0, 0, 0;
     L2(2,1), L2(2,2) + L3(1,1), L3(1,2), 0;
     0, L3(2,1), L3(2,2), 0;
     L1(2,1), 0, 0, L1(2,2)];
 
 
 Yinv = inv(Y);
 
 
 for i = 1:5000
 
    napon_novi = Yinv * I;
    
    delta = 0.1;
    
    p_sen = p;
    p_sen(4) = p_sen(4) + j*delta;
    
    I_sen = izracunajStruje(p_sen, napon);
    napon_sen = Yinv * I_sen;
    
    sensitivity = (abs(napon_sen(4)) - abs(napon(4))) / delta;
    
    promijena = j*(abs(napon_novi(4)) - abs(napon(4))) / sensitivity;
    
    p(4) = p(4) + promijena;

    I = izracunajStruje(p, napon_novi);
    
    if abs(promijena) < 1e-6
        napon = napon_novi;
        napon_novi = Yinv * I;
        disp("interacija: ");
        disp(i);
        break;
    end
    
    junk = napon;
    napon = napon_novi;

 end
 disp("napon: ");
 disp(napon);
 disp("provjera ispravnosti (Y*napon - I (ako su naponi ispravni treba da bude jednako 0))");
 disp(Y*napon - I);

 
 