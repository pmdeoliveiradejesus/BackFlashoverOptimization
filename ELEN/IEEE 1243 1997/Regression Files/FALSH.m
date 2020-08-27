%function [BFR] = FLASH(Vnom,X,Y,topcrossarm,TR,R,Dbundle,span,W,FR,GFD)
%% IEEE 1243 - Guide for Improving the Lightning Performance of Transmission Lines
%% Input parameters 
    Vnom = 345; %kV
    % coordinates;
    X = [-5.5 5.5 -5.5 -8.6  -5.8 5.5 8.6  5.8]; %m
    Y = [39.3 39.3 33.8 27.4  21.3 33.8 27.4 21.3];  %m
    R = [0.45 0.45 1.48 1.48 1.48 1.48 1.48 1.48 ]/100;  % m radius
    Dbundle = 0.467;    %m distance in the bundle
    topcrossarm = [2.7 9.3 15.3];   %[m]
    span = 335; % [m] 
    W = 2.63;   % [m] Insulation Length 
    GFD = 3.6;  % Flash density per km2 per year
    FR = 20;    % Flooting resistance
    TR = 10;    % [m] Tower Base Radius
%% IEEE 1243 - Guide for Improving the Lightning Performance of Transmission Lines
    %% Setup Data
    %% Coordinates 
    % [m]  Conductor Horizontal Distance
    Xg = X(1:2);    % Shield Wire
    X = X(3:8);
    % [m] Conductor Vertical Distance
    Yg = Y(1:2);   % Shield Wire
    Y = Y(3:8);  % [m]
    TW = Yg(1);
    BW = Xg(2)-Xg(1); % Horizontal separation between two shield wires
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
    %GFD = 3.6;
    RN = GFD*(28*TW^(0.6)+BW)/10; %Line Flashes per 100km per year
    % Effective tower flashes
    EF = 0.6 * RN;
    %% Tower Modeling 
    TS = SP / (VL * 0.9);

    % 2.1 and 2.2 - flashover voltages at 2 and 6 us
    if TS >= 10
        F2 = W * (400 + 710 * (2 * TS) ^ -0.75);
    else
        F2 = 820 * W;
    end
    F6 = 585 * W;

    % 2.3 - tower-top voltage
    TV = 1.8 * F2; % TODO - always the first phase conductor?

    % 2.4 - shield wire corona adjustment
    RC = fsolve(@(RC) RC*log((2*TW)/RC)-TV/Eo,[1],optimoptions('fsolve','Display','off')); % TODO - always the first shield wire?

    % 2.5 - self surge impedance of each shield wire
    GZ = 60 * sqrt(log(2 * TW / R(1)) * log(2 * TW / RC));

    % 2.6 - combined surge impedance
    if NE == 1
        GC = GZ;
    else  % TODO - hard wired for 2 shield wires
        BM = abs(Xg(2)-Xg(1));
        AM = sqrt((2*TW)^2+BM^2);
        GM = 60 * log(AM / BM);
        GC = (GZ + GM) / 2;
    end

    % 2.7 - mutual impedances between conductors and shield wires
    for i = 1 : NE
        for j = 1 : NC
            AM = sqrt((TW+Y(j))^2+(Xg(i)-X(j))^2);
            BM = sqrt((TW-Y(j))^2+(Xg(i)-X(j))^2);
            CZ(j, i) = 60 * log(AM / BM);
        end
    end

    % 2.7b - coupling factors
    for i = 1 : NC
        if NE == 1
            CF(i) = CZ(i, 1) / GC;
        else   % TODO - hard wired for 2
            CF(i) = (CZ(i, 1) + CZ(i, 2)) / 2 / GC;
        end
    end

    % 2.8 and 2.9 - tower surge response
    TU = TW / (VL * 0.85);
    ZT = 30 * log(2 * (TW^2 + TR * TR / 4) / (TR * TR / 4));

    % 2.11 - travel time to each crossarm
    for i = 1 : NC
        TA(i) = (TW - Y(i)) / (VL * 0.85);  % TODO - parameterize tower wave velocity
    end

    % 2.12 and 2.13 - intrinsic circuit impedance
    ZI = (GC * ZT) / (GC + 2 * ZT);

    % start the main footing resistance loop
    %FR = 20;

    % 2.14 - tower wave impedance
    ZW = (2 * ZT * GC * GC) / (GC + 2 * ZT)^2 * (ZT - FR) / (ZT + FR);

    % 2.15 - tower damping factor
    PS = (2 * ZT - GC) / (2 * ZT + GC) * (ZT - FR) / (ZT + FR);
    % betaS is not used
    % betaS = (2 * ZT - GC) / (2 * ZT + GC)


    % 2.16 - footing resistance damping factor
    AP = (2 * FR) / (ZT + FR);

    % 2.17 through 2.19 - tower top voltage
    V2 = ZI - ZW * (1 - TU / (1 - PS)) / (1 - PS);
    if TS < 1
        VT = -4 * 0.85 * V2 * V2 / GC * (1 - 2 * V2 / GC) * (1 - TS);
    else
        VT = 0;
    end
    V2 = V2 + VT;

    % 2.20 - voltage across footing resistance at 2 us
    VW = AP * ZI * (1 - (PS * TU) / (1 - PS)) / (1 - PS);

    % 2.21 - reduce the voltage across footing resistance by reflection
    VW = VW + VW * VT / V2;

    % 2.22 - crossarm voltages at 2 us
    for i = 1 : NC
        VC(i) = VW + ((TU - TA(i)) / TU) * (V2 - VW);
    end

    % 2.23 - insulator voltages at 2 us
    for i = 1 : NC
        VI(i) = VC(i) - CF(i) * V2;
    end

    % 2.24 - tower top voltage at 6 us
    V6 = GC * FR / (GC + 2 * FR);

    % 2.25 - reflected voltage at 6 us
    D9 = FR / (GC + 2 * FR);
    VV = -4 * 0.85 * GC * D9 * D9 * (1 - 2 * D9);

    % 2.26 - insulator voltages at 6 us
    for i = 1 : NC
        VS(i) = (V6 + VV) * (1 - CF(i));
    end

    % 2.27 - critical stroke current at 2 us
    for i = 1 : NC
        I2(i) = F2 / VI(i);
    end

    % 2.28 - critical stroke current at 6 us
    for i = 1 : NC
        I6(i) = F6 / VS(i);

        % 2.29 and 2.30 - select lowest critical current
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
    
TH = [0 -2*pi/3 2*pi/3 2*pi/3 -2*pi/3 0];   % theta

for i = 1:NC
    IC2(i) = (VR(i) - Vnom*cos(TH(i))) / VI(i);
    P(i) = 1/(1+(IC2(i)/31)^2.6);
end
disp('NL')
disp(RN)
PBF = max(P);
disp('Probability')
disp(sum(PBF)/length(PBF))
disp('BFR')
BFR = EF*sum(PBF)/length(PBF)