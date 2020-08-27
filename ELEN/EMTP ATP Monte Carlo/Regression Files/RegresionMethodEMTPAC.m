clear
clc
close all
%% Datos
BFR = [];
R = 5:5:50;
W = 200:25:600;
u = 1:1:1;
GFD=3.6;
BW = 11;
TW = 39.3;
RN = 0.1 * GFD * (BW + 28*TW^0.6);
BFR=0.6*RN*[
0.0043200000000000	0.0021800000000000	0.0013600000000000	0.0008200000000000	0.0003600000000000	0.0002200000000000	0.0001200000000000	0.0000600000000000	0.0000400000000000	0.0000200000000000	0.0000200000000000	0.0000200000000000	0.0000200000000000	0.0000200000000000	0.0000133333333333	0.0000066666666667	0.0000033666666667	...
0.0146428571428571	0.0077571428571429	0.0045571428571429	0.0028571428571429	0.0017142857142857	0.0011142857142857	0.0007285714285714	0.0005000000000000	0.0003571428571429	0.0002285714285714	0.0001285714285714	0.0001142857142857	0.0000428571428571	0.0000428571428571	0.0000285714285714	0.0000142857142857	0.0000102857142857	...
0.0367000000000000	0.0197000000000000	0.0127000000000000	0.0075000000000000	0.0053000000000000	0.0039000000000000	0.0026000000000000	0.0014000000000000	0.0014000000000000	0.0009000000000000	0.0006000000000000	0.0004000000000000	0.0002000000000000	0.0001500000000000	0.0001000000000000	0.0000800000000000	0.0000500000000000	...
0.0700000000000000	0.0422000000000000	0.0281000000000000	0.0169000000000000	0.0125000000000000	0.0077000000000000	0.0056000000000000	0.0043000000000000	0.0034000000000000	0.0022000000000000	0.0015000000000000	0.0013000000000000	0.0009000000000000	0.0006000000000000	0.0005000000000000	0.0004000000000000	0.0002000000000000	...
0.1039000000000000	0.0700000000000000	0.0462000000000000	0.0327000000000000	0.0211000000000000	0.0146000000000000	0.0099000000000000	0.0070000000000000	0.0054000000000000	0.0045000000000000	0.0036000000000000	0.0023000000000000	0.0015000000000000	0.0013000000000000	0.0009000000000000	0.0009000000000000	0.0005000000000000	...
0.1443000000000000	0.0990000000000000	0.0696000000000000	0.0475000000000000	0.0345000000000000	0.0247000000000000	0.0171000000000000	0.0135000000000000	0.0096000000000000	0.0075000000000000	0.0056000000000000	0.0045000000000000	0.0035000000000000	0.0030000000000000	0.0025000000000000	0.0024000000000000	0.0014000000000000	...
0.1915000000000000	0.1351000000000000	0.0941000000000000	0.0698000000000000	0.0496000000000000	0.0374000000000000	0.0281000000000000	0.0203000000000000	0.0145000000000000	0.0117000000000000	0.0088000000000000	0.0070000000000000	0.0053000000000000	0.0044000000000000	0.0037000000000000	0.0030000000000000	0.0028000000000000	...
0.2075000000000000	0.1468333333333330	0.1091666666666670	0.0820000000000000	0.0590000000000000	0.0431666666666667	0.0305000000000000	0.0251666666666667	0.0183333333333333	0.0133333333333333	0.0103333333333333	0.0080000000000000	0.0060000000000000	0.0043333333333333	0.0038333333333333	0.0033333333333333	0.0033333333333333	...
0.2913057961359090	0.2210193204530310	0.1670552964690210	0.1240839440373080	0.0911059293804131	0.0699533644237175	0.0532978014656895	0.0396402398401066	0.0321452365089940	0.0251499000666223	0.0191538974017322	0.0143237841439041	0.0119920053297801	0.0098267821452365	0.0074950033311126	0.0059960026648901	0.0048301132578281	...
0.3426000000000000	0.2659000000000000	0.2028000000000000	0.1528000000000000	0.1159000000000000	0.0925000000000000	0.0700000000000000	0.0538000000000000	0.0408000000000000	0.0333000000000000	0.0267000000000000	0.0209000000000000	0.0160000000000000	0.0135000000000000	0.0101000000000000	0.0087000000000000	0.0073000000000000	]'/GFD;     




%% Variables de ajuste por m?nimos cuadrados
c0 = ones(length(R)*length(W),1);
c1 = repmat(R,length(W),1);
c1 = c1(:);
c2 = W.*ones(1,length(W)); 
c2 = repmat(c2(:),1,length(R)).*ones(length(W),length(R)); 
c2 = c2(:);
c3 = ones(1,length(W)*length(R)); 
c3 = c3(:);
cx3=repmat(c3,1,7);
cx2=repmat(c2,1,7);
cx1=repmat(c1,1,7);
 X = [[c0 c1 c1.^2 c1.^3 c1.^4 c1.^5 c1.^6].*cx3 [c0 c1 c1.^2 c1.^3 c1.^4 c1.^5 c1.^6].*cx2.*cx3 ...
     [c0 c1 c1.^2 c1.^3 c1.^4 c1.^5 c1.^6].*cx2.^2.*cx3 [c0 c1 c1.^2 c1.^3 c1.^4 c1.^5 c1.^6].*cx2.^3.*cx3...
     [c0 c1 c1.^2 c1.^3 c1.^4 c1.^5 c1.^6].*cx2.^4.*cx3 [c0 c1 c1.^2 c1.^3 c1.^4 c1.^5 c1.^6].*cx2.^5.*cx3...
     [c0 c1 c1.^2 c1.^3 c1.^4 c1.^5 c1.^6].*cx2.^6.*cx3];
%% Ajuste por mInimos cuadrados
theta = (X'*X)\X'*BFR;
R2 = 1 - sum((BFR - X*theta).^2)/sum((BFR - mean(BFR)).^2); 
disp(R2);
save thetaEMTPAC.mat theta
Rg=5;
We=500;
BFRest =GFD*[[1 Rg Rg^2 Rg^3 Rg^4 Rg^5 Rg^6] ...
     [1 Rg Rg^2 Rg^3 Rg^4 Rg^5 Rg^6]*(We) ...
     [1 Rg Rg^2 Rg^3 Rg^4 Rg^5 Rg^6]*(We)^2 ...
     [1 Rg Rg^2 Rg^3 Rg^4 Rg^5 Rg^6]*(We)^3 ...
     [1 Rg Rg^2 Rg^3 Rg^4 Rg^5 Rg^6]*(We)^4 ...
     [1 Rg Rg^2 Rg^3 Rg^4 Rg^5 Rg^6]*(We)^5 ...
     [1 Rg Rg^2 Rg^3 Rg^4 Rg^5 Rg^6]*(We)^6]*theta
 



