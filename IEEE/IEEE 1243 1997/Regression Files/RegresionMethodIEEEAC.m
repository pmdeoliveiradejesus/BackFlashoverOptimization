clear
clc
close all

%% Datos
 BFR = [];
R = 1:1:50;
W = 200:1:600;
u = 1:1:1;
for i = 1:length(R)
    for j = 1:length(W)
        for k = 1:length(u)
            BFR = [BFR ; IEEE_method(R(i),W(j),u(k))]; %M?todo de c?lculo de tasa de salida
        end
    end
    i
end
save BFR BFR
%load('BFR')
%% Variables de ajuste por m?nimos cuadrados
% c0 = ones(length(R)*length(W)*length(u),1);
% c1 = R.*ones(length(W)*length(u),length(R)); c1 = c1(:);
% c2 = W.*ones(length(u),length(W)); c2 = c2(:).*ones(length(u)*length(W),length(R)); c2 = c2(:);
% c3 = (u').*ones(length(u),length(W)*length(R)); c3 = c3(:);
% X = [[c0 c1 c1.^2 c1.^3 c1.^4 c1.^5 c1.^6].*c3 [c0 c1 c1.^2 c1.^3 c1.^4 c1.^5 c1.^6].*c2.*c3 ...
%     [c0 c1 c1.^2 c1.^3 c1.^4 c1.^5 c1.^6].*c2.^2.*c3 [c0 c1 c1.^2 c1.^3 c1.^4 c1.^5 c1.^6].*c2.^3.*c3...
%     [c0 c1 c1.^2 c1.^3 c1.^4 c1.^5 c1.^6].*c2.^4.*c3 [c0 c1 c1.^2 c1.^3 c1.^4 c1.^5 c1.^6].*c2.^5.*c3...
%     [c0 c1 c1.^2 c1.^3 c1.^4 c1.^5 c1.^6].*c2.^6.*c3];
% 
% %% Ajuste por m?nimos cuadrados
% theta = (X'*X)\X'*BFR;
% R2 = 1 - sum((BFR - X*theta).^2)/sum((BFR - mean(BFR)).^2); 
% disp(R2);
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
save theta theta
% % %% Graficar BFR vs Rg con W y u constante
%   %  Rg = 10:1:20;
%     BFOR = [];
%     x = [];
%     i = 20;
%     k = 3.6;
% for j = 1:length(W)
%     BFOR = [BFOR ; IEEE_method(i,W(j),k)];
%     x = [x ; IEEE_methodEst(i,W(j),k)];
% end
% figure
% plot(W,BFOR);
% title('BFOR vs W');
% grid on 
% hold on 
% plot(W,x ); 
% %% Graficar BFR vs Rg con W y u constante
% %Rg = 10:1:20;
%     BFOR = [];
%     y = [];
%     j = 263;
%     k = 3.6;
% for i = 1:length(R)
%     BFOR = [BFOR ; IEEE_method(R(i),j,k)];
%     y = [y ; IEEE_methodEst(R(i),j,k)];
% end
% figure
% plot(R,BFOR);
% title('BFOR vs Rg');
%  grid on 
%  hold on 
%  plot(R,y );

