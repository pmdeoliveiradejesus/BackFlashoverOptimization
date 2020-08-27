clear
clc
close all
%% Datos
BFR = [];
R = 1:1:50;
W = 200:2:600;
u = 1:1:1;
for i = 1:length(R)
    for j = 1:length(W)
        for k = 1:length(u)
            BFR = [BFR ; EpriBFRAC(R(i),W(j),u(k))]; %M?todo de c?lculo de tasa de salida
        end
    end
    disp(i);
end
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
%% Ajuste por m?nimos cuadrados
theta = (X'*X)\X'*BFR;
R2 = 1 - sum((BFR - X*theta).^2)/sum((BFR - mean(BFR)).^2); 
disp(R2);
save thetaEpriAC.mat theta
load('thetaEpriAC')%load fit
GFD=3.6;
Rg=20;
W=275;
BFRest =GFD*[[1 Rg Rg^2 Rg^3 Rg^4 Rg^5 Rg^6] ...
     [1 Rg Rg^2 Rg^3 Rg^4 Rg^5 Rg^6]*(W) ...
     [1 Rg Rg^2 Rg^3 Rg^4 Rg^5 Rg^6]*(W)^2 ...
     [1 Rg Rg^2 Rg^3 Rg^4 Rg^5 Rg^6]*(W)^3 ...
     [1 Rg Rg^2 Rg^3 Rg^4 Rg^5 Rg^6]*(W)^4 ...
     [1 Rg Rg^2 Rg^3 Rg^4 Rg^5 Rg^6]*(W)^5 ...
     [1 Rg Rg^2 Rg^3 Rg^4 Rg^5 Rg^6]*(W)^6]*theta