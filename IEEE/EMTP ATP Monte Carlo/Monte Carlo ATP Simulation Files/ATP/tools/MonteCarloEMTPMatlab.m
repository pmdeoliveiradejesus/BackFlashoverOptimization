%function [BFR] = ATP(Rg,w,GFD)
%%
clc
clear all
%pkg load statistics
%% Datos de diseño
    delete('BFR.pl4');
    delete('BFR.lis');
    delete('BFR.atp'); 
    delete('bfr01.mat');
Rg=1;%Resistencia de Impulso de Puesta a Tierra, ohms 
GFD=3.6;% Densidad de rayos, rayos/km^2/año
w=[200 225 250 275 300 325 350 375 400 425 450 475 500 525 550 575 600];% Aislamiento en cm
%w=263;
D = 1; %Cantidad de descargas aleatorias;
%%% Otros datos
BW = 11; % distancia entre cables de guarda
TW = 39.3; % altura de la torre
%RN = 0.1 * GFD * (BW + 28*TW^0.6); %frecuencia de impacto, IEEE 1247
RN = 0.1 * GFD * (BW + 4*TW^1.09); %frecuencia de impacto, EPRI Red Book
%Inicializacion del contador de descargas exitosas
   for  h=1:length(w) 
nfo(1,h) = 0; 
   end
% determinación de la curva CFO-tiempo en KV   
t=1:1:8001;% %nanosegundos
for  h=1:length(w)   
     CFO(:,h) = 0.01*w(h) * (400 + 710 * (t/1000) .^ -0.75) * 1000; %Calcula el CFO in V
    end
kk=1;
for i = 1:D
    % angulo de la fuente de tensión aleatorio en grados
    ang=round(unifrnd(1,360));
    angA=(ang+0);
    angB=(ang-120);    
    angC=(ang-240);
    I =  round(lognrnd(log(33300),0.605)); %Corriente aleatoria de la descarga, kA       
    tf = round(lognrnd(log(2),0.494),1)*10^-6; %Tiempo de frente de onda de la descarga, microsegundos
    %Modificación del archivo de ATP para incluir la nueva corriente,
    %tiempo de frente de onda de la descarga y Rg.
 fid = fopen('BFRx','w');
    fprintf(fid,horzcat('BEGIN NEW DATA CASE\n$DUMMY, XYZ000\n   1.E-9   8.E-6                        \n     500       1       1       1       1       0       0       1       0\n/BRANCH\n-1XX0029XX0007                    027.2.55E8   2.7 1 0                         0\n-1XX0007XX0006                    117.2.55E8   6.6 1 0                         0\n-1XX0006XX0005                    111.2.55E8    6. 1 0                         0\n-1XX0005XX0004                    116.2.55E8   24. 1 0                         0\n  XX0004                     ',num2str(Rg),'.                                               0\n-1XX0030XX0011                    027.2.55E8   2.7 1 0                         0\n-1XX0011XX0010                    117.2.55E8   6.6 1 0                         0\n-1XX0010XX0009                    111.2.55E8    6. 1 0                         0\n-1XX0009XX0008                    116.2.55E8   24. 1 0                         0\n  XX0008                     ',num2str(Rg),'.                                               0\n-1XX0031PhaseA                    027.2.55E8   2.7 1 0                         0\n-1PhaseAPhaseB                    117.2.55E8   6.6 1 0                         0\n-1PhaseBPhaseC                    111.2.55E8    6. 1 0                         0\n-1PhaseCXX0012                    116.2.55E8   24. 1 0                         0\n  XX0012                     ',num2str(Rg),'.                                               0\n-1XX0032XX0016                    027.2.55E8   2.7 1 0                         0\n-1XX0016XX0015                    117.2.55E8   6.6 1 0                         0\n-1XX0015XX0014                    111.2.55E8    6. 1 0                         0\n-1XX0014XX0013                    116.2.55E8   24. 1 0                         0\n  XX0013                     ',num2str(Rg),'.                                               0\n-1XX0033XX0020                    027.2.55E8   2.7 1 0                         0\n-1XX0020XX0019                    117.2.55E8   6.6 1 0                         0\n-1XX0019XX0018                    111.2.55E8    6. 1 0                         0\n-1XX0018XX0017                    116.2.55E8   24. 1 0                         0\n  XX0017                     ',num2str(Rg),'.                                               0\n$INCLUDE, C:/ATP/work/span.lib, PHA###, PHB###, PHC###, PH2A##, PH2B## $$\n  , PH2C##, XX0031, XX0031, X0026A, X0026B, X0026C, X0025A, X0025B, X0025C $$\n  , XX0032, XX0032\n$INCLUDE, C:/ATP/work/span.lib, X0024A, X0024B, X0024C, X0023A, X0023B $$\n  , X0023C, XX0030, XX0030, PHA###, PHB###, PHC###, PH2A##, PH2B##, PH2C## $$\n  , XX0031, XX0031\n$INCLUDE, C:/ATP/work/span.lib, X0026A, X0026B, X0026C, X0025A, X0025B $$\n  , X0025C, XX0032, XX0032, X0028A, X0028B, X0028C, X0027A, X0027B, X0027C $$\n  , XX0033, XX0033\n$INCLUDE, C:/ATP/work/spanf.lib, X0028A, X0028B, X0028C, X0027A, X0027B $$\n  , X0027C, XX0033, XX0033, X0001A, X0001B, X0001C, X0001A, X0001B, X0001C $$\n  , ######, ######\n$INCLUDE, C:/ATP/work/span.lib, X0022A, X0022B, X0022C, X0021A, X0021B $$\n  , X0021C, XX0029, XX0029, X0024A, X0024B, X0024C, X0023A, X0023B, X0023C $$\n  , XX0030, XX0030\n$INCLUDE, C:/ATP/work/spanf.lib, X0002A, X0002B, X0002C, X0003A, X0003B $$\n  , X0003C, ######, ######, X0022A, X0022B, X0022C, X0021A, X0021B, X0021C $$\n  , XX0029, XX0029\n/SOURCE\n12XX0031-1 ',           num2str(I,'%6.2e'),'              ',num2str(tf,'%1.1e'),'                                  1.\n14X0001A   281691.32       60. ',num2str(angA,'%5.2e'),'                        -1       100.\n14X0001B   281691.32       60. ',num2str(angB,'%5.2e'),'                       -1       100.\n14X0001C   281691.32       60. ',num2str(angC,'%5.2e'),'                       -1       100.\n/OUTPUT\n  PHA   PHB   PHC   PhaseAPhaseBPhaseCPH2A  PH2B  PH2C  \nBLANK BRANCH\nBLANK SWITCH\nBLANK SOURCE\nBLANK OUTPUT\nBLANK PLOT\nBEGIN NEW DATA CASE\nBLANKK CASE\nBLANK'));   
  fclose(fid); 
    fid = fopen('BFR.ATP','w');
    
    %fprintf(fid,horzcat('BEGIN NEW DATA CASE\n$DUMMY, XYZ000\n   1.E-9   8.E-6                        \n     500       1       1       1       1       0       0       1       0\n/BRANCH\n-1XX0029XX0007                    145.2.55E8   2.7 1 0                         0\n-1XX0007XX0006                    145.2.55E8   6.6 1 0                         0\n-1XX0006XX0005                    145.2.55E8    6. 1 0                         0\n-1XX0005XX0004                    145.2.55E8   24. 1 0                         0\n  XX0004                     ',num2str(Rg),'.                                               0\n-1XX0030XX0011                    145.2.55E8   2.7 1 0                         0\n-1XX0011XX0010                    145.2.55E8   6.6 1 0                         0\n-1XX0010XX0009                    145.2.55E8    6. 1 0                         0\n-1XX0009XX0008                    145.2.55E8   24. 1 0                         0\n  XX0008                     ',num2str(Rg),'.                                               0\n-1XX0031PhaseA                    145.2.55E8   2.7 1 0                         0\n-1PhaseAPhaseB                    145.2.55E8   6.6 1 0                         0\n-1PhaseBPhaseC                    145.2.55E8    6. 1 0                         0\n-1PhaseCXX0012                    145.2.55E8   24. 1 0                         0\n  XX0012                     ',num2str(Rg),'.                                               0\n-1XX0032XX0016                    145.2.55E8   2.7 1 0                         0\n-1XX0016XX0015                    145.2.55E8   6.6 1 0                         0\n-1XX0015XX0014                    145.2.55E8    6. 1 0                         0\n-1XX0014XX0013                    145.2.55E8   24. 1 0                         0\n  XX0013                     ',num2str(Rg),'.                                               0\n-1XX0033XX0020                    145.2.55E8   2.7 1 0                         0\n-1XX0020XX0019                    145.2.55E8   6.6 1 0                         0\n-1XX0019XX0018                    145.2.55E8    6. 1 0                         0\n-1XX0018XX0017                    145.2.55E8   24. 1 0                         0\n  XX0017                     ',num2str(Rg),'.                                               0\n$INCLUDE, C:/ATP/work/span.lib, PHA###, PHB###, PHC###, PH2A##, PH2B## $$\n  , PH2C##, XX0031, XX0031, X0026A, X0026B, X0026C, X0025A, X0025B, X0025C $$\n  , XX0032, XX0032\n$INCLUDE, C:/ATP/work/span.lib, X0024A, X0024B, X0024C, X0023A, X0023B $$\n  , X0023C, XX0030, XX0030, PHA###, PHB###, PHC###, PH2A##, PH2B##, PH2C## $$\n  , XX0031, XX0031\n$INCLUDE, C:/ATP/work/span.lib, X0026A, X0026B, X0026C, X0025A, X0025B $$\n  , X0025C, XX0032, XX0032, X0028A, X0028B, X0028C, X0027A, X0027B, X0027C $$\n  , XX0033, XX0033\n$INCLUDE, C:/ATP/work/spanf.lib, X0028A, X0028B, X0028C, X0027A, X0027B $$\n  , X0027C, XX0033, XX0033, X0001A, X0001B, X0001C, X0001A, X0001B, X0001C $$\n  , ######, ######\n$INCLUDE, C:/ATP/work/span.lib, X0022A, X0022B, X0022C, X0021A, X0021B $$\n  , X0021C, XX0029, XX0029, X0024A, X0024B, X0024C, X0023A, X0023B, X0023C $$\n  , XX0030, XX0030\n$INCLUDE, C:/ATP/work/spanf.lib, X0002A, X0002B, X0002C, X0003A, X0003B $$\n  , X0003C, ######, ######, X0022A, X0022B, X0022C, X0021A, X0021B, X0021C $$\n  , XX0029, XX0029\n/SOURCE\n12XX0031-1 ',           num2str(I,'%6.2e'),'              ',num2str(tf,'%1.1e'),'                                  1.\n14X0001A   281691.32       60. ',num2str(angA,'%5.2e'),'                        -1       100.\n14X0001B   281691.32       60. ',num2str(angB,'%5.2e'),'                       -1       100.\n14X0001C   281691.32       60. ',num2str(angC,'%5.2e'),'                       -1       100.\n/OUTPUT\n  PHA   PHB   PHC   PhaseAPhaseBPhaseCPH2A  PH2B  PH2C  \nBLANK BRANCH\nBLANK SWITCH\nBLANK SOURCE\nBLANK OUTPUT\nBLANK PLOT\nBEGIN NEW DATA CASE\nBLANKK CASE\nBLANK'));   
  
    fprintf(fid,horzcat('BEGIN NEW DATA CASE\n$DUMMY, XYZ000\n   1.E-9   8.E-6                        \n     500       1       1       1       1       0       0       1       0\n/BRANCH\n-1XX0029XX0007                    027.2.55E8   2.7 1 0                         0\n-1XX0007XX0006                    117.2.55E8   6.6 1 0                         0\n-1XX0006XX0005                    111.2.55E8    6. 1 0                         0\n-1XX0005XX0004                    116.2.55E8   24. 1 0                         0\n  XX0004                     ',num2str(Rg),'.                                               0\n-1XX0030XX0011                    027.2.55E8   2.7 1 0                         0\n-1XX0011XX0010                    117.2.55E8   6.6 1 0                         0\n-1XX0010XX0009                    111.2.55E8    6. 1 0                         0\n-1XX0009XX0008                    116.2.55E8   24. 1 0                         0\n  XX0008                     ',num2str(Rg),'.                                               0\n-1XX0031PhaseA                    027.2.55E8   2.7 1 0                         0\n-1PhaseAPhaseB                    117.2.55E8   6.6 1 0                         0\n-1PhaseBPhaseC                    111.2.55E8    6. 1 0                         0\n-1PhaseCXX0012                    116.2.55E8   24. 1 0                         0\n  XX0012                     ',num2str(Rg),'.                                               0\n-1XX0032XX0016                    027.2.55E8   2.7 1 0                         0\n-1XX0016XX0015                    117.2.55E8   6.6 1 0                         0\n-1XX0015XX0014                    111.2.55E8    6. 1 0                         0\n-1XX0014XX0013                    116.2.55E8   24. 1 0                         0\n  XX0013                     ',num2str(Rg),'.                                               0\n-1XX0033XX0020                    027.2.55E8   2.7 1 0                         0\n-1XX0020XX0019                    117.2.55E8   6.6 1 0                         0\n-1XX0019XX0018                    111.2.55E8    6. 1 0                         0\n-1XX0018XX0017                    116.2.55E8   24. 1 0                         0\n  XX0017                     ',num2str(Rg),'.                                               0\n$INCLUDE, C:/ATP/work/span.lib, PHA###, PHB###, PHC###, PH2A##, PH2B## $$\n  , PH2C##, XX0031, XX0031, X0026A, X0026B, X0026C, X0025A, X0025B, X0025C $$\n  , XX0032, XX0032\n$INCLUDE, C:/ATP/work/span.lib, X0024A, X0024B, X0024C, X0023A, X0023B $$\n  , X0023C, XX0030, XX0030, PHA###, PHB###, PHC###, PH2A##, PH2B##, PH2C## $$\n  , XX0031, XX0031\n$INCLUDE, C:/ATP/work/span.lib, X0026A, X0026B, X0026C, X0025A, X0025B $$\n  , X0025C, XX0032, XX0032, X0028A, X0028B, X0028C, X0027A, X0027B, X0027C $$\n  , XX0033, XX0033\n$INCLUDE, C:/ATP/work/spanf.lib, X0028A, X0028B, X0028C, X0027A, X0027B $$\n  , X0027C, XX0033, XX0033, X0001A, X0001B, X0001C, X0001A, X0001B, X0001C $$\n  , ######, ######\n$INCLUDE, C:/ATP/work/span.lib, X0022A, X0022B, X0022C, X0021A, X0021B $$\n  , X0021C, XX0029, XX0029, X0024A, X0024B, X0024C, X0023A, X0023B, X0023C $$\n  , XX0030, XX0030\n$INCLUDE, C:/ATP/work/spanf.lib, X0002A, X0002B, X0002C, X0003A, X0003B $$\n  , X0003C, ######, ######, X0022A, X0022B, X0022C, X0021A, X0021B, X0021C $$\n  , XX0029, XX0029\n/SOURCE\n12XX0031-1 ',           num2str(I,'%6.2e'),'              ',num2str(tf,'%1.1e'),'                                  1.\n14X0001A   281691.32       60. ',num2str(angA,'%5.2e'),'                        -1       100.\n14X0001B   281691.32       60. ',num2str(angB,'%5.2e'),'                       -1       100.\n14X0001C   281691.32       60. ',num2str(angC,'%5.2e'),'                       -1       100.\n/OUTPUT\n  PHA   PHB   PHC   PhaseAPhaseBPhaseCPH2A  PH2B  PH2C  \nBLANK BRANCH\nBLANK SWITCH\nBLANK SOURCE\nBLANK OUTPUT\nBLANK PLOT\nBEGIN NEW DATA CASE\nBLANKK CASE\nBLANK'));   
    fclose(fid);
    system ('runATP.vbs'); %Corre el archivo de ATP a partir de un archivo .vbs  
    system('pl42mat.vbs'); %Convierte los resultados de ATP de .pmatlab alll4 a .mat a través de un .vbs y GTPPL
    load('bfr01.mat');  
    %Carga el archivo .mat
        VI = vPHASEA - vPHA; 
    VI = [VI vPHASEB - vPHB];
    VI = [VI vPHASEC - vPHC];
    VI = [VI vPHASEC - vPH2A];
    VI = [VI vPHASEB - vPH2B];
    VI = [VI vPHASEA - vPH2C];
        for  h=1:length(w)   
    %Calcula la tension en la cadena de aisladores para cada fase
    %Verifica si hay backflashover en cada fase
    IC(:,1) = VI(:,1) > CFO(:,h);
    IC(:,2) = VI(:,2) > CFO(:,h);
    IC(:,3) = VI(:,3) > CFO(:,h);
    IC(:,4) = VI(:,4) > CFO(:,h);
    IC(:,5) = VI(:,5) > CFO(:,h);
    IC(:,6) = VI(:,6) > CFO(:,h);
    nfo(kk+1,h) = nfo(kk,h) + (sum(sum(IC)) > 0);   
    P(kk,h) = nfo(kk+1,h)/kk; %Probabilidad de ocurrencia de backflashover
T(kk,h) = 0.6*RN*P(kk,h);
    end

[   
%     w(1)/100   kk nfo(kk,1) P(kk,1) T(kk,1);
%      w(2)/100   kk nfo(kk,2) P(kk,2) T(kk,2);
%      w(3)/100   kk nfo(kk,3) P(kk,3) T(kk,3);
     w(4)/100   kk/1000 nfo(kk,4) P(kk,4) T(kk,4);
%      w(5)/100   kk nfo(kk,5) P(kk,5) T(kk,5);
]

kk=kk+1;
end









