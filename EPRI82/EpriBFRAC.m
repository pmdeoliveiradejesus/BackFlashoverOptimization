function [Tbf] = EpriBFRAC(Rg,W,Ng)
%% work sheet 1-A 345kV Vertical double circuit two ground wires
%% Setup tower coordinates and additional data
W=W/100;%string length in meters
Vnom = 345; %kV
% coordinates;
X = [-5.5 5.5 -5.5 -8.6  -5.8 5.5 8.6  5.8]; %m
Y = [39.3 39.3 33.8 27.4  21.3 33.8 27.4 21.3];  %m
R = [0.45 0.45 1.48 1.48 1.48 1.48 1.48 1.48 ]/100;  % m radius
Dbundle = 0.467;    %m distance in the bundle
%Rg = 20;    % [Ohm]
%W = 2.63; % m Insulatorlength
span = 335; %m 
topcrossarm = [2.7 9.3 15.3];   %[m]
TR = 10;    % [m] Tower Base Radius
T = 30; % Kereunic Level
Vo = sqrt(2)*Vnom/sqrt(3);  %kV
    hs = 7;   %m midspan
    Eo = 1500;    % kV/m 
    top2crosam1 = topcrossarm(1);    %m
    top2crosam2 = topcrossarm(2);    %m
    top2crosam3 = topcrossarm(3);   %m
    Po = 760; %mmHg, atmospheric pressure at sea level
    H = 0;    %m, Altitude
    Ta = 20;  %Celsius degrees, air temperature
    C = 0.003671; %Thermal expansion coefficient 1/C
    P = Po/10^(H/(18400*(1+C*Ta))); %Presure mmHg
    Kt = 0.3855*P/(273+Ta);   % Correction pressure and temperature
    %% Step 2 Define or determine Keraunic level
    %% Step 3 Define or determine Stroke incidence to earth
    %N = 0.12*T; %flash density per km2 per year Eq. 12.4.1
    %% Step 4 Compute mean shield wire height
    hg = Y(1)-(2/3)*(hs); %m, equivalent height of shield wires Eq. 12.4.3
    %% Step 5 Compute total flashes to the line NL=Nsf+Nbf
    NL = Ng*((X(2)-X(1))+4*hg^1.09)/10;  % Line Flashes per 100km per year Eq. 12.4.6 (Epri)
    %NL = Ng*(28*Y(1)^(0.6)+X(2)-X(1))/10;  % Line Flashes per 100km per year (IEEE1243) 
    %% Step 6 Look up Flashover Voltage of insulation strings at each phase
    % Equations in Fig 12.6.3
    t = 6;  %mus
    K1 = 0.4*W;
    K2 = 0.71*W;
    Vs6 = 1000*(K1+K2/t^0.75)*Kt;   %kV V impulse Volt-time curve 
    %% Step 7 Compute mean phase wire height
    for j = 3:8
        hf(j) = Y(j)-(2/3)*hs; %m Eq. 12.4.3
    end
    %% Step 8 Compute single radius in all phases with Corona effect
    for j=3:8
        x = fsolve(@(x) x*log((2*hf(j))/x)-Vs6/Eo,[1],optimoptions('fsolve','Display','off'));
        rFcor(j) = x; %m Radius of corona envelope Eq. 12.5.2
    end
    %% Step 9 Compute equivalent single radius in all phases without Corona effect
    for j=3:8
        req(j) = sqrt(R(j)*Dbundle); %m Equivalent radius Eq. 12.5.1
    end
    %% Step 10 Compute correted radius in all phases without Corona effect
    for j=3:8
        RC(j) = req(j)+rFcor(j);%m
    end
    %% Step 11 Compute effective self-surge impedance in all phase
    for j = 3:8
        Zf(j) = 60*sqrt(log(2*hf(j)/req(j))*log(2*hf(j)/RC(j))); %ohms Eq. 12.5.3
    end
    %% Step 12 Compute minimum stroke current for shielding failure Eq. 12.7.6
    for j = 3:8
        Imin(j) = 2*Vs6/Zf(j);  %kA
    end
    %% Step 13 Compute minimum distance for shielding failure
    for j=3:8
        Smin(j) = 10*Imin(j)^0.65;    %m Eq. 12.7.1
    end
    %% Step 14 Select beta value
    beta = 0.8;
    % Step 15 Omitted
    % Step 16 Omitted
    %% Computation of actual shielding failure rate 
    %% Step 17 Compute uncovered width Xs per phase (<0 perfect shielded) 
    % Eq. 12.7.2 || Eq. 12.7.3 
    Yg = hg;
    for j=3:5
        Yf(j) = hf(j);
        Xf(j) = X(j);
        Xg = X(1);

        Xp(j) = sqrt(Smin(j)^2-(beta*Smin(j)-Yf(j))^2)+Xf(j);
        Xge(j) = -sqrt(Smin(j)^2-(beta*Smin(j)-Yg)^2)+Xp(j);
        F(j) = sqrt((Yg-Yf(j))^2+(Xg-Xf(j))^2);
        theta(j) = asin((beta*Smin(j)-Yf(j))/Smin(j));
        w(j) = acos(F(j)/(2*Smin(j)));
        alphas(j) = abs(atan((Xf(j)-Xg)/(Yg-Yf(j))));

        if Yf(j)-beta*Smin(j) > 0
            theta(j)=1;
        end
        Xs(j) = subplus(Smin(j)*(cos(theta(j))+sin(alphas(j)-w(j))));% m
    end

    for j=6:8
        Yf(j) = hf(j);
        Xf(j) = X(j);
        Xg = X(2);
        Xp(j) = sqrt(Smin(j)^2-(beta*Smin(j)-Yf(j))^2)+Xf(j);
        Xge(j) = -sqrt(Smin(j)^2-(beta*Smin(j)-Yg)^2)+Xp(j);
        F(j) = sqrt((Yg-Yf(j))^2+(Xg-Xf(j))^2);
        theta(j) = asin((beta*Smin(j)-Yf(j))/Smin(j));
        w(j) = acos(F(j)/(2*Smin(j)));
        alphas(j) = atan((Xf(j)-Xg)/(Yg-Yf(j)));

        if Yf(j)-beta*Smin(j)>0
            theta(j) = 1;
        end
        Xs(j) = Smin(j)*(cos(theta(j))+sin(alphas(j)-w(j))); %m
    end

    Xs(Xs<0) = 0; % only positive widths
    %% Step 18 Compute maximun strike distance that can occur per phase, Eq. 12.7.7
    for j = 3:5
        Xg = X(1);
        Yo(j) = (Yf(j)+Yg)/2;
        m(j) = abs((Xf(j)-Xg)/(Yg-Yf(j)));
        As(j) = m(j)^2-beta*m(j)-beta^2;
        Bs(j) = beta*(m(j)^2+1);
        Cs(j) = (m(j)^2+1);
        Sp(j) = (-Bs(j)-sqrt(Bs(j)^2-4*As(j)*Cs(j)))/(2*As(j));
        Smax(j) = Yo(j)*Sp(j);
    end
    for j = 6:8
        Xg = X(2);
        Yo(j) = (Yf(j)+Yg)/2;
        m(j) = abs((Xf(j)-Xg)/(Yg-Yf(j)));
        As(j) = m(j)^2-beta*m(j)-beta^2;
        Bs(j) = beta*(m(j)^2+1);
        Cs(j) = (m(j)^2+1);
        Sp(j) = (-Bs(j)-sqrt(Bs(j)^2-4*As(j)*Cs(j)))/(2*As(j));
        Smax(j) = Yo(j)*Sp(j);% m
    end
    %Step 19 ---> Step 12 
    %% Step 20 Compute Imax Eq. 12.7.1A
    for j = 3:8
        Imax(j) = 0.029*Smax(j)^1.54;% kA
    end
    %% Step 21 and 22 FIND PROBABILITY THAT Imin and Imax will be exceeded
    for j = 3:8
        ProbMINsf(j) = 1/(1+(Imin(j)/31)^2.6);
        ProbMAXsf(j) = 1/(1+(Imax(j)/31)^2.6);
    end
    %% Step 23 Compute shielding failures per circuit and phase
    for j = 3:8
        Nsf(j) = (Ng/10)*(Xs(j))*(ProbMINsf(j)-ProbMAXsf(j))/2; %shielding failures/100km*yr
    end
    %% Step 24 Compute total shielding failures
    NSft = sum(Nsf);  %shielding failures/100km/yr %outages/100km.yr (shield failure)
    %% Step 25 
    Nbf = NL-NSft; %Total flashes/100km/yr to be used in backflashover assessmen
    %% -------------------------------------------------------------------------------- %%
    %% ------- Back-flashover calculations (EPRI 1982) Two Point Method --------------- %%
    %% -------------------------------------------------------------------------------- %%
    %% STEP 1: Insulator flasover voltage kV at 2 museg
    t = 2;  %museg
    K1 = 0.4*W;
    K2 = 0.71*W;
    Vs2 = 1000*(K1+K2/t^0.75)*Kt;   %KV
    Vs2x = 820*W*Kt; %Eq. 12.10.10
    %% STEP 2: Insulator flasover voltage kV at 6museg
    Vs6;    % Step 6 
    Vs6x = 585*W*Kt;  %Eq. 12.10.11
    %% STEP 3: Multiply Vs2 by 1.8
    Vsn2 = 1.8*Vs2; % Tower top voltage and average for all phase 
    %% STEP 4: COMPUTE SHIELD WIRE CORONA RADIUS
    % Corona effect in shield wire
    x = fsolve(@(x) x*log((2*Y(1))/x)-Vsn2/Eo,[1],optimoptions('fsolve','Display','off'));
    Rcor = x;%m Figure 12.5.3 Eo=1500kv/m Eq. 12.5.2
    %% STEP 5 Compute self-surge and mutial impedance of each shield wire
    Zcor(1,1) = 60*sqrt(log(2*Y(1)/R(1))*log(2*Y(1)/Rcor)); %with high corona Eq. 12.5.3
    %% Step 6 Compute combined surge impedance (two shield wires)
    b(1,2) = abs(X(1)-X(2));
    a(1,2) = sqrt(b(1,2)^2+(2*Y(1))^2);
    %Z(1,1) = 60*log(2*Y(1)/R(1));% ohms
    Z(1,2) = 60*log(a(1,2)/b(1,2)); % ohms Eq. 12.5.5A
    Zs = (Z(1,2)+Zcor(1,1))/2;% ohm  Eq. 12.5.6
    %% Step 7 Compute coupling factors 
    % Calculate self and mutual surge impedances of all phase, ohm
    for j=3:8
        b(j,1) = sqrt((X(1)-X(j))^2+(Y(1)-Y(j))^2);
        b(j,2) = sqrt((X(2)-X(j))^2+(Y(2)-Y(j))^2);
        a(j,1) = sqrt(b(j,1)^2+(Y(j)+Y(1))^2);
        a(j,2) = sqrt(b(j,2)^2+(Y(j)+Y(2))^2);
        Z(j,1) = 60*log(a(j,1)/b(j,1)); % Eq 12.5.5A
        Z(j,2) = 60*log(a(j,2)/b(j,2));
        k(j) = (Z(j,1)+Z(j,2))/(Z(1,2)+Zcor(1,1));  % coupling factor Eq. 12.5.7
    end
    %% Step 8 Compute Tower Surge Impedance
    % Figure 12.5.5
    h = Y(1);
    r = 10/2;
    Zt = 30*log(2*(h^2+r^2)/r^2);
    %% Step 9 Determine tower travel time
    taot = Y(1)/300;    % mus 
    %% Step 10 Determine span travel time
    taoS = span/(300*0.9);  % mus
    %% Step 11 Travel time tpn from tower top to each crossarm  
    tao(3) = (top2crosam1)/(300);   % mus 
    tao(4) = (top2crosam2)/(300);   % mus 
    tao(5) = (top2crosam3)/(300);   % mus 
    tao(6) = (top2crosam1)/(300);   % mus 
    tao(7) = (top2crosam2)/(300);   % mus 
    tao(8) = (top2crosam3)/(300);   % mus 
    %% STEP 12 set impulse footing resistance Rg
    %ohm impulse footing resistance Fig. 12.6.1
    %% STEP 13 intrinsic circuit impedance
    ZI = Zs*Zt/(Zt*2+Zs);   %ohm   Eq. 12.6.1A
    %% STEP 14 tower wave impedance
    Zw =((2*Zs^2*Zt)/(Zt*2+Zs)^2)*((Zt-Rg)/(Zt+Rg));    %ohm   Eq. 12.6.1B
    %% STEP 15 tower damping
    psi =((2*Zt-Zs)/(Zt*2+Zs))*((Zt-Rg)/(Zt+Rg));%    Eq. 12.6.1C
    %% STEP 16 footing reflection factor
    alphaR =((2*Rg)/(Zt+Rg));%    Eq. 12.6.4
    %% STEP 17 per unit tower top voltage at 2 mus
    Vt2 = (ZI-(Zw/(1-psi))*(1-taot/(1-psi)));%    Eq. 12.10.1
    %% STEP 18 reflected component from adjacent towers at 2mus
    Ks = 0.85;
    if taoS > 1
        Vtp2 = 0;%  NO reflection  Eq. 12.10.3
    else
        Vtp2 = (-4*Ks*Vt2^2/Zs)*((1-2*Vt2)/Zs)*(1-taoS); % Eq. 12.10.3
    end
    %% STEP 19 actual tower top voltage at 2mus
    Vtt2 = Vt2+Vtp2;
    %% STEP 20 voltage across footing resistance at 2mus
    Vr2 = (alphaR*ZI/(1-psi))*(1-psi*taot/(1-psi));%    Eq. 12.10.2
    %% STEP 21 reflected component from adjacent towers at 2mus
    Ks = 0.85;
    if taoS > 1
        Vrp2 = 0;%  NO reflection  Eq. 12.10.3
    else
        Vrp2 = (-4*Ks*Vr2^2/Zs)*((1-2*Vr2)/Zs)*(1-taoS) ;%    Eq. 12.10.3
    end
    Vtr2 = Vr2+Vrp2;
    %% STEP 22 compute crossarm voltage at 2mus (at each phase) Eq. 12.10.5
    for j=3:8
        Vn2(j) = Vtr2+((taot-tao(j))/taot)*(Vtt2-Vtr2);
    end
    %% STEP 23 compute per unit insulator voltage at 2mus (at each phase) Eq. 12.10.6
    for j=3:8
        Vsn2(j) = Vn2(j)-k(j)*Vtt2;
    end
    %% STEP 24 compute tower top voltage at 6mus  Eq. 12.10.7
    Vt6 = (Zs*Rg)/(Zs+2*Rg);%  
    %% STEP 25 compute tower top voltage at 6mus (reflected)   
    Vtp6 = (-4*Ks*Zs)*(Rg/(Zs+2*Rg))^2*(1-2*Rg/(Zs+2*Rg)) ;%    Eq. 12.10.8
    %% STEP 26 compute per unit insulator voltage at 6mus (at each phase) Eq. 12.10.9
    for j=3:8
        Vsn6(j) = (Vt6+Vtp6)*(1-k(j));
    end
    %% STEP 27 compute ratios (Icn)2 at 2mus (at each phase) in kA Eq. 12.10.12
    for j = 3:8
        Icn2(j) = (Vs2)/Vsn2(j);
    end
    %% STEP 28 compute ratios (Icn)6 at 6mus (at each phase) in kA Eq. 12.10.13
    for j=3:8
        Icn6(j) = (Vs6)/Vsn6(j);
    end
    %% Step 29 select the lowest currents (all at 2mus)
    %% Step 30 select the corresponding Vs 
    Vs2;
    %% Step 31 plot IC fluctuation, percentage of time
    n = 360;
    for theta=1:n
        Icreal(3,theta) = Icn2(3)*(Vs2-Vo*sin(pi*theta/180-0))/Vs2;   %Phase A
        Icreal(4,theta) = Icn2(4)*(Vs2-Vo*sin(pi*theta/180+2*pi/3))/Vs2;  %Phase B
        Icreal(5,theta) = Icn2(5)*(Vs2-Vo*sin(pi*theta/180-2*pi/3))/Vs2;  %Phase C
        Icreal(6,theta) = Icn2(6)*(Vs2-Vo*sin(pi*theta/180-2*pi/3))/Vs2;  %Phase C'
        Icreal(7,theta) = Icn2(7)*(Vs2-Vo*sin(pi*theta/180+2*pi/3))/Vs2;  %Phase B'
        Icreal(8,theta) = Icn2(8)*(Vs2-Vo*sin(pi*theta/180+0))/Vs2;   %Phase A'
        th(theta) = theta;
    end
    %% Step 32 
    for theta = 1:n
        A(theta) = min([Icreal(3,theta) Icreal(4,theta) Icreal(5,theta) Icreal(6,theta) Icreal(7,theta) Icreal(8,theta) ]);
        for j = 3:8
           if A(theta) == Icreal(j,theta)
                B(theta) = j;
           end
        end
    end

    for j=3:8
        ki(j) = 0;
    end

    for theta = 1:n
        for j = 3:8
            if B(theta)== j
                ki(j)= ki(j)+1;
            end  
        end
    end

    for j=3:8
        p(j)=ki(j)/n;
    end

    for j=3:8
        Ix(j)=0;
    end

    for theta=1:n
        for j=3:8
            if B(theta) == j
                Ix(j) = Ix(j)+A(theta);
            end
        end 
    end
    %% Step 33 Critical flashover current kA
    for j=3:8
        Ic(j)=Ix(j)/(ki(j)+.0001);
    end
    %% Step 34
    for j=3:8
        ProbBF(j)=1/(1+(Ic(j)/31)^2.6);% Probability that stroke current Ic will be exceeded in any flash to the line
    end
    %% Step 35 Effective tower flashes per 100 km per year 
    eff = 0.6; % 60% 
    NE = eff*(Nbf);% Effective tower flashes per 100km per yr
    %% Step 36 Tower flashes per phase 
    Np = NE*p; %flashes/100km.yr (per phase)
    %% Step 37 Expected  number of strokes causing flashover of a given phase 
    Tbfp = Np.*ProbBF;%outages/100km.yr (only backflshover, per phase)
    %% Step 38 Total backflashover per 100 km per year 
% %     disp('Total backflashover per 100 km per year')
% %     disp('RED Book')
    Tbf = sum(Tbfp); %outages/100km.yr (only backflshover)
    % %% step 39 Total failures per 100 km per year 
    % disp('Total failures per 100 km per year')
    % Tbfplussh = Tbf+NSft %outages/100km.yr (backflshover and shield failure)
%     %% Figures 
%     figure
%     plot(th,Icreal(3,:),'r','LineWidth',1)
%     hold on 
%     plot(th,Icreal(4,:),'b','LineWidth',1) 
%     plot(th,Icreal(5,:),'g','LineWidth',1) 
%     plot(th,Icreal(6,:),'g--','LineWidth',1)
%     plot(th,Icreal(7,:),'b--','LineWidth',1)
%     plot(th,Icreal(8,:),'r--','LineWidth',1)
% 
%     title('IC fluctuation per each phase - EPRI')
%     ylabel('Icn [kA]')
%     xlabel('angle [?]')
%     legend('Phase A','Phase B','Phase C','Phase C"','Phase B"','Phase A"')
%     hold off
% 
%     figure
%     plot(th,A,'r','LineWidth',1)
%     title('Average Icn per phase dominated - EPRI')
%     xlabel('angle [?]')
%     ylabel('Icn [kA]')
%     axis([0 360 122 136])
%     xticks(0:30:360)
% 
% 
%     figure
%     plot(th,B-3*ones(1,length(B)),'r','LineWidth',1)
%     title('Phase Dominated - EPRI')
%     xlabel('angle [?]')
%     ylabel('phase')
%     axis([0 360 1 6])
%     xticks(0:30:360)
%     yticks(1:1:6)
%     yticklabels({'A','B','C','C"','B"','A"'})
end

