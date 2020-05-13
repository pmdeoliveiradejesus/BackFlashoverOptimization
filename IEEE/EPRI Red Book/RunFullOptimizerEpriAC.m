% This program solves otlgimulse vaying Tmax from 0.1 to 1.075
% in order to build Figure 11 
clear all
close all
clc
global theta b kI rho Tmax density CG CI L n Lt Tx w60  wmax R0   NPH NT bw Eog Irayo Rb NL
load('thetaEpriAC'); %Load fit model
b=[	55.52532416	 	1.221891523	];
kI=25;%incremental tower insulation cost $/cm
BW = 11;
TW = 39.3;
bw = (BW + 28*TW^0.6);
Eog=500;%kV/m
Lt=100;%km
w60=263;%60Hz and switching insulation
wmax=445;%cm
R0=50000000;%ohms default tower resistance
n=12;%number of sections
rho=[	850	750	1250	1150  975        1125        1500        2400	 4000	4200	2000	1800	];
density=[	4.6	4.6	4.6	4.6		3.6	3.6	3.6	3.6	2.6	2.6	2.6	2.6 ];
NPH=ones(1,n)*6;% Number of phases per section
NT=ones(1,n)*300/n; %Number of towers per section
for k=1:n
L(k)=Lt/n;%section length km
end																								
x0a=19.9626*ones(1,n);
x0b=263.0000*ones(1,n);
tol1=1e-10;  
tol2=1e-10;
tol3=1e-10;
  LB=[];
  UB=[];
A=[];
B=[];
Aeq=[];
Beq=[];
for u=1:1:160
Tmax=.1+u/100;
%Tmax=1.4;%BF tripout rate 
x0=horzcat(x0a,x0b);
%options=optimset('Display','iter','LargeScale','on','ActiveConstrTol',1);
options=optimset('Algorithm','interior-point','Display','iter','TolFun',tol1,'TolCon',tol2,'TolX',tol3,'MaxIter',30000000,'MaxFunEvals',10000000);
%FMINCON Call objetive function and non-linear constraints 
[x,fval,exitflag,output,lambda,grad,hessian]=fmincon('minCostimp6fit',x0,A,B,Aeq,Beq,LB,UB,'restrCost2imp6fit',options);
stat(u)=exitflag;
%% Determine TOTAL AND Incremental costs
mu=lambda.ineqnonlin';
lmd=lambda.eqnonlin;
for k=1:n%  
Tx(k) =density(k)*[[1 x(k) x(k)^2 x(k)^3 x(k)^4 x(k)^5 x(k)^6] ...
     [1 x(k) x(k)^2 x(k)^3 x(k)^4 x(k)^5 x(k)^6]*(x(k+n)) ...
     [1 x(k) x(k)^2 x(k)^3 x(k)^4 x(k)^5 x(k)^6]*(x(k+n))^2 ...
     [1 x(k) x(k)^2 x(k)^3 x(k)^4 x(k)^5 x(k)^6]*(x(k+n))^3 ...
     [1 x(k) x(k)^2 x(k)^3 x(k)^4 x(k)^5 x(k)^6]*(x(k+n))^4 ...
     [1 x(k) x(k)^2 x(k)^3 x(k)^4 x(k)^5 x(k)^6]*(x(k+n))^5 ...
     [1 x(k) x(k)^2 x(k)^3 x(k)^4 x(k)^5 x(k)^6]*(x(k+n))^6]*theta;
 NL(k)=density(k)*(bw)/10;
Irayo(k)=31*((0.6*NL(k))/Tx(k)-1)^(1/2.6);
Rb(k)=sqrt((rho(k)*Eog*x(k)^2)/(rho(k)*Eog-2*pi*Irayo(k)*x(k)^2));   
CG(k)=NT(k)*b(1)*(rho(k)/Rb(k))^(b(2));
CI(k)=NT(k)*NPH(k)*kI*(x(k+n)-w60);
CG1(u)=sum(CG)+0*sum(CI);
CI1(u)=0*sum(CG)+sum(CI);
Ctot(u)=sum(CG)+sum(CI);
Cgtot(u)=+sum(CG);
syslambda(u)=-lmd;
countl(u)=u;
if stat(u) < 0
Ctot(u)=Ctot(u-1);
syslambda(u)=syslambda(u-1);
end
end
end
close all
figure
plot(Ctot)
figure
plot(syslambda)
