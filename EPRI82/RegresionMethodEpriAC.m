clear
clc
close all
%% Datos
BFR = [];
R = 1:.1:50;
W = 200:1:600;
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
c0 = ones(length(R)*length(W)*length(u),1);
c1 = R.*ones(length(W)*length(u),length(R)); c1 = c1(:);
c2 = W.*ones(length(u),length(W)); c2 = c2(:).*ones(length(u)*length(W),length(R)); c2 = c2(:);
c3 = (u').*ones(length(u),length(W)*length(R)); c3 = c3(:);
X = [[c0 c1 c1.^2 c1.^3 c1.^4 c1.^5 c1.^6].*c3 [c0 c1 c1.^2 c1.^3 c1.^4 c1.^5 c1.^6].*c2.*c3 ...
    [c0 c1 c1.^2 c1.^3 c1.^4 c1.^5 c1.^6].*c2.^2.*c3 [c0 c1 c1.^2 c1.^3 c1.^4 c1.^5 c1.^6].*c2.^3.*c3...
    [c0 c1 c1.^2 c1.^3 c1.^4 c1.^5 c1.^6].*c2.^4.*c3 [c0 c1 c1.^2 c1.^3 c1.^4 c1.^5 c1.^6].*c2.^5.*c3...
    [c0 c1 c1.^2 c1.^3 c1.^4 c1.^5 c1.^6].*c2.^6.*c3];
%% Ajuste por m?nimos cuadrados
theta = (X'*X)\X'*BFR;
R2 = 1 - sum((BFR - X*theta).^2)/sum((BFR - mean(BFR)).^2); 
disp(R2);
% %% Graficar BFR vs Rg con W y u constante
% Rg = 1:50;
% BFOR = [];
% x = [];
% j = 2.63;
% k = 3.6;
% for i = 1:length(Rg)
%     BFOR = [BFOR ; EpriRedAC2(Rg(i),j,k)];
%     x = [x; [1 i i^2 i^3 i^4 i^5 i^6]*k [1 i i^2 i^3 i^4 i^5 i^6]*j*k [1 i i^2 i^3 i^4 i^5 i^6]*j^2*k ...
%             [1 i i^2 i^3 i^4 i^5 i^6]*j^3*k [1 i i^2 i^3 i^4 i^5 i^6]*j^4*k [1 i i^2 i^3 i^4 i^5 i^6]*j^5*k...
%             [1 i i^2 i^3 i^4 i^5 i^6]*j^6*k];
% end
% plot(Rg,BFOR);
% title('BFOR vs Rg');
% grid on
% hold on
% plot(Rg,x*theta);
% %% Graficar BFR vs W con Rg y u constante
% W = 2:0.1:6;
% BFOR = [];
% x = [];
% i = 20;
% k = 3.6;
% for j = 2:0.1:6
%     BFOR = [BFOR ; CIGRE2(i,j,k)];
%     x = [x; [1 i i^2 i^3 i^4 i^5 i^6]*k [1 i i^2 i^3 i^4 i^5 i^6]*j*k [1 i i^2 i^3 i^4 i^5 i^6]*j^2*k ...
%             [1 i i^2 i^3 i^4 i^5 i^6]*j^3*k [1 i i^2 i^3 i^4 i^5 i^6]*j^4*k [1 i i^2 i^3 i^4 i^5 i^6]*j^5*k...
%             [1 i i^2 i^3 i^4 i^5 i^6]*j^6*k];
% end
% plot(W,BFOR);
% title('BFOR vs W');
% grid on
% hold on
% plot(W,x*theta);
% %% Graficar BFR vs u con Rg y W constante
% W = 2:0.1:12;
% BFOR = [];
% x = [];
% i = 20;
% j = 2.63;
% for k = 2:0.1:12
%     BFOR = [BFOR ; EpriRedAC(i,j,k)];
%     x = [x; [1 i i^2 i^3 i^4 i^5 i^6]*k [1 i i^2 i^3 i^4 i^5 i^6]*j*k [1 i i^2 i^3 i^4 i^5 i^6]*j^2*k ...
%             [1 i i^2 i^3 i^4 i^5 i^6]*j^3*k [1 i i^2 i^3 i^4 i^5 i^6]*j^4*k [1 i i^2 i^3 i^4 i^5 i^6]*j^5*k...
%             [1 i i^2 i^3 i^4 i^5 i^6]*j^6*k];
% end
% plot(W,BFOR);
% title('BFOR vs \mu');
% grid on
% hold on
% plot(W,x*theta);
save theta theta
