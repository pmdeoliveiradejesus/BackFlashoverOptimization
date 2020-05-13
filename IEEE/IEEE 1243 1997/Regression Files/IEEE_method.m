function [BFR] = IEEE(FR,W,GFD)
    %% IEEE 1243 - Guide for Improving the Lightning Performance of Transmission Lines
    %% Setup Data
Vnom = 345; %kV
% coordinates;
W=W/100;
X = [-5.5 5.5 -5.5 -8.6  -5.8 5.5 8.6  5.8]; %m
Y = [39.3 39.3 33.8 27.4  21.3 33.8 27.4 21.3];  %m
R = [0.45 0.45 1.48 1.48 1.48 1.48 1.48 1.48 ]/100;  % m radius
Dbundle = 0.467;    %m distance in the bundle
span = 335; %m 
topcrossarm = [2.7 9.3 15.3];   %[m]
TR = 10;    % [m] Tower Base Radius
%T = 30; % Kereunic Level

    
    %% Coordinates 
    % [m]  Conductor Horizontal Distance
    Xg = X(1:2);    % Shield Wire
    X = X(3:8);
    % [m] Conductor Vertical Distance
    Yg = Y(1:2);   % Shield Wire
    Y = Y(3:8);  % [m]
    TW = Yg(1);
    BW = Xg(2)-Xg(1); % Horizontal separation between two shield wires
    % Input Values 
    % TR --> Tower Base Radius
    % FR --> Flooting Resistance
    % T --> keraunic level tunderstorm days per year
    % W --> [m] Insulator Length
    %% Tower Data 
    SP = span;   % Span [m]
    NC = 6; % Number of conductors
    NE = 2; % Number of shield wire 
    NB = 2; % number of bundled subconductors
    BS = Dbundle; %m BundleSpacing

    Eo = 1500;  % kV/m 
    VL = 300;   % m/mus
    t2 = 2;  % muS
    t6 = 6;  % muS

    top2crosam1 = topcrossarm(1);  % [m]
    top2crosam2 = topcrossarm(2);  % [m]
    top2crosam3 = topcrossarm(3); % [m]
    MAX_ANGLES = 180;
    MAX_PHASES = 6;
    PA = [0 240 120 120 240 0];
    PV = sqrt(2)/sqrt(3);
    LV = Vnom*ones(1,6);    % Line Voltage
    %% Lightning 
    %GFD = 0.04*T^1.25;   % Ground Flash Density
    RN = GFD*(28*TW^(0.6)+BW)/10; %Line Flashes per 100km per year

    % Effective tower flashes
    EF = 0.6 * RN;
    %% Tower Modeling 
    % span travel time
    TS = SP / (VL * 0.9);
    % Flashover voltages at 2 and 6 us
    if TS >= 10
        F2 = W * (400 + 710 * (2 * TS) ^ -0.75);
    else
    F2 = 820 * W;
    end
    F6 = 585 * W;

    % Tower-top voltage
    TV = 1.8 * F2;
    % Shield wire corona adjustment
    RC = fsolve(@(RC) RC*log((2*TW)/RC)-TV/Eo,[1],optimoptions('fsolve','Display','off')); 

    % Self surge impedance of each shield wire
    GZ = 60 * sqrt(log(2 * TW / R(1)) * log(2 * TW / RC));

    % Combined surge impedance
    if NE == 1
        GC = GZ;
    else
        BM = abs(Xg(2)-Xg(1));
        AM = sqrt((2*TW)^2+BM^2);
        GM = 60 * log(AM / BM);
        GC = (GZ + GM) / 2;
    end

    % Mutual impedances between conductors and shield wires
    for i = 1 : NE
        for j = 1 : NC
            AM = sqrt((TW+Y(j))^2+(Xg(i)-X(j))^2);
            BM = sqrt((TW-Y(j))^2+(Xg(i)-X(j))^2);
            CZ(j, i) = 60 * log(AM / BM);
        end
    end

    % Coupling factors
    for i = 1 : NC
        if NE == 1
            CF(i) = CZ(i, 1) / GC;
        else   
            CF(i) = (CZ(i, 1) + CZ(i, 2)) / 2 / GC;
        end
    end

    % Tower surge response
    TU = TW / (VL * 0.85);
    ZT = 30 * log(2 * (TW^2 + TR * TR / 4) / (TR * TR / 4));

    % Travel time to each crossarm
    for i = 1 : NC
        TA(i) = (TW - Y(i)) / (VL * 0.85);
    end
    % Intrinsic circuit impedance
    ZI = (GC * ZT) / (GC + 2 * ZT);

    % Start the main footing resistance loop
    % Footing Resistance

    % Tower wave impedance
    ZW = (2 * ZT * GC * GC) / (GC + 2 * ZT)^2 * (ZT - FR) / (ZT + FR);

    % Tower damping factor
    PS = (2* ZT - GC) / (2 * ZT + GC) * (ZT - FR) / (ZT + FR);
    % betaS = (2 * ZT - GC) / (2 * ZT + GC)

    % Footing resistance damping factor
    AP = (2 * FR) / (ZT + FR);

    % Tower top voltage
    V2 = ZI - ZW * (1 - TU / (1 - PS)) / (1 - PS);
    if TS < 1
        VT = -4 * 0.85 * V2 * V2 / GC * (1 - 2 * V2 / GC) * (1 - TS);
    else
        VT = 0;
    end
    V2 = V2 + VT;

    % Voltage across footing resistance at 2 us
    VW = AP * ZI * (1 - (PS * TU) / (1 - PS)) / (1 - PS);
    % Reduce the voltage across footing resistance by reflection
    VW = VW + VW * VT / V2;

    % Crossarm voltages at 2 us
    for i = 1 : NC
        VC(i) = VW + ((TU - TA(i)) / TU) * (V2 - VW);
    end

    % 2.23 - insulator voltages at 2 us
    for i = 1 : NC
        VI(i) = VC(i) - CF(i) * V2;
    end

    % Tower top voltage at 6 us
    V6 = GC * FR / (GC + 2 * FR);

    % Reflected voltage at 6 us
    D9 = FR / (GC + 2 * FR);
    VV = -4 * 0.85 * GC * D9 * D9 * (1 - 2 * D9);

    % Insulator voltages at 6 us
    for i = 1 : NC
        VS(i) = (V6 + VV) * (1 - CF(i));
    end
    %% Critical Stroke Current 
    % Critical stroke current at 2 us
    for i = 1 : NC
        I2(i) = F2 / VI(i);
    end
    % Critical stroke current at 6 us
    for i = 1 : NC
        I6(i) = F6 / VS(i);
        % Select lowest critical current
        if I2(i) < I6(i)
            ISCRIT(i) = I2(i);
            VR(i) = F2;
        else
            ISCRIT(i) = I6(i);
            VR(i) = F6;
        end
        IU(i) = round(ISCRIT(i),2);
        VJ(i) = round(VI(i),2);
        VP(i) = round(VS(i),2);
    end

    %% Find equation for critical current
    d_ang = 360 / MAX_ANGLES;
    ang = 0.5 * d_ang;
    for i = 1 : MAX_ANGLES
        SV(i) = sin(ang * pi / 180);
        ang = ang + d_ang;
    end
    for i = 1 : MAX_PHASES
        AD(i) = 0;
        DN(i) = 0;
        PR(i) = 0;
        for j = 1 : MAX_ANGLES
            KS(i, j) = 0;
            KE(i, j) = 0;
        end
    end
    d_k = 360 / MAX_ANGLES;
    k = d_k / 2;
    while k < 360
        IL = 2000; %initial guess of min current (kA)
        for j = 1 : NC
            CC = (k + PA(j)) / d_k;
            CC = int32(CC);
            while CC < 1
                CC = CC + MAX_ANGLES;
            end
            while CC > MAX_ANGLES
                CC = CC - MAX_ANGLES;
            end
            PH(j) = SV(CC);
            IR(j) = ISCRIT(j) * (1 - PV * PH(j) * LV(j) / VR(j));
            if IR(j) <= IL
                IL = IR(j);
            end
        end
        %Flash 1.9 - all phases with IR <= IL share dominance equally
        NUM_IL = 0;
        for j = 1 : NC
            if IR(j) <= IL
                NUM_IL = NUM_IL + 1;
            end
        end
        if NUM_IL > 0
            share_IL = 1 / NUM_IL;
            i_ang = k - d_k / 2;
            for j = 1 : NC
                if IR(j) <= IL % start or extend a dominance interval
                    PR(j) = PR(j) + share_IL;
                    if DN(j) < 1  % start a new interval
                        B = AD(j) + 1;
                        KS(j, B) = i_ang;
                        DN(j) = 1;
                        AD(j) = AD(j) + 1;
                    end
                    B = AD(j);
                    KE(j, B) = i_ang + d_k;
                else   %end a dominance interval
                    DN(j) = 0;
                end
            end
        end
        k = k + d_k;
    end

    % Percent of time each phase dominates
    for i = 1 : NC
        PR(i) = PR(i) * 100 / MAX_ANGLES;
    end

    % Average value of current for each phase
    for i = 1 : NC
        SU(i) = 0;
        SA(i) = 0;
        jmax = AD(i);  %the number of dominant intervals to average
        if jmax < 1
            IC(i) = 0;
        end
        if jmax >=1
            for j = 1 : jmax
                T2 = (KE(i, j) + PA(i)) * pi / 180;
                T1 = (KS(i, j) + PA(i)) * pi / 180;
                TH(j) = T2 - T1;
                IA(j) = ISCRIT(i) * (1 + (PV * LV(i) / VR(i)) * (cos(T2) - cos(T1)) / TH(j));
                SU(i) = SU(i) + IA(j) * TH(j);
                SA(i) = SA(i) + TH(j);
            end
            IC(i) = SU(i) / SA(i);
        end
    end

    % Probability that stroke current is exceeded
    for i = 1 : NC
        if IC(i) > 0
            PX(i) = 1/(1+(IC(i)/31)^2.6);
        else
            PX(i) = 0;
        end
    end
   % 2.36 through 2.38 - tower flashes per phase per 100 km per year
    SS = 0;
    for i = 1 : NC
        ET(i) = EF * PR(i) * PX(i) / 100;
        SS = SS + ET(i);
    end
    SS = round(SS,2);

    TH = [0 -2*pi/3 2*pi/3 2*pi/3 -2*pi/3 0];   % theta
    for i = 1:NC
        IC2(i) = (VR(i) - 345*cos(TH(i))) / VI(i);
        P(i) = 1/(1+(IC2(i)/31)^2.6);
    end
    Icritical = min(IC2);
    PBF = max(P);
    BFR = EF*sum(PBF)/length(PBF);
%     disp('Total backflashover per 100 km per year - IEEE1243')
%     disp( VI)
end

