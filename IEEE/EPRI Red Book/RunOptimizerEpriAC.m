% NEW OTLGI
% incuding CIGRE resistance impulse formula
clear all
close all
clc
global theta b kI rho Tmax density CG CI L n Lt Tx w60  wmax R0 NPH NT bw Eog Irayo Rb NL
load('thetaEpriAC')%load fit
b=[	55.52532416	 	1.221891523	];%Groundig cost parameters
kI=25;%incremental tower insulation cost $/cm
bw1=201.5938;%m
BW = 11;
TW = 39.3;
bw = (BW + 28*TW^0.6);
Eog=500;%kV/m
Lt=100;%km
Tmax=1.4;%BF tripout rate 
w60=263;%60Hz and switching insulation
wmax=445;%cm
R0=50000000;%ohms

% % Case 1 One section
% n=1;%number of sections
%  rho=[1833];
%  density=[3.6];

% 
% %Case 2 3 sections
% rho=[1000 1500 3000];
% density=[4.6 3.6 2.6];

%Case 3 12 towers
n=12;
rho=[	850	750	1250	1150  975        1125        1500        2400	 4000	4200	2000	1800	];
density=[	4.6	4.6	4.6	4.6		3.6	3.6	3.6	3.6	2.6	2.6	2.6	2.6 ];

% 
% 
% %Case 4 300 towers
%  rhox=[	850	750	1250	1150  975        1125        1500        2400	 4000	4200	2000	1800	];
%  densityx=[	4.6	4.6	4.6	4.6		3.6	3.6	3.6	3.6	2.6	2.6	2.6	2.6 ];
% j=0;
% e=.10;
% e2=0.02;
% for i=1:8  
% for k=1:25
%     j=j+1;
%     rho(j)=(1-e/2)*rhox(i)+(k-1)*(e)*rhox(i)/24;      %0.01*rhox(i)*(-1 + (2)*rand(1,1));
%     density(j)=densityx(i);%+0.01*densityx(i)*(-1 + (2)*rand(1,1));
% end
% end
% for i=9:12  
% for k=1:25
%     j=j+1;
%     rho(j)=(1.0+e2/2)*rhox(i)-(k-1)*(e2)*rhox(i)/24;      %0.01*rhox(i)*(-1 + (2)*rand(1,1));
%     density(j)=densityx(i);%+0.01*densityx(i)*(-1 + (2)*rand(1,1));
% end
% end
% %close all
% %plot(rho)
% e3=0;
% j=0;
% for i=1:12  
% for k=1:25
%     j=j+1;
%     density(j)=(1-e3/2)*densityx(i)+(k-1)*(e3)*densityx(i)/24;      %0.01*rhox(i)*(-1 + (2)*rand(1,1));
%     %density(j)=densityx(i);%+0.01*densityx(i)*(-1 + (2)*rand(1,1));
% end
% end




NPH=ones(1,n)*6;% Number of phases per section
NT=ones(1,n)*300/n; %Number of towers per section
for k=1:n
L(k)=Lt/n;%section length km
end																								
  

x0a=19.9626*ones(1,n);
x0b=263.0000*ones(1,n);
%x0a=30*ones(1,n);
%x0b=400.0000*ones(1,n);
tol1=1e-6;  
tol2=1e-6;
tol3=1e-6;
  LB=[];
  UB=[];
A=[];
B=[];
Aeq=[];
Beq=[];
x0=horzcat(x0a,x0b);
%for u=1:1;  
%Tmax=u/10;
%options=optimset('Display','iter','LargeScale','on','ActiveConstrTol',1);
options=optimset('Algorithm','interior-point','Display','iter','TolFun',tol1,'TolCon',tol2,'TolX',tol3,'MaxIter',30000000,'MaxFunEvals',10000000);
%FMINCON Call objetive function and non-linear constraints 
[x,fval,exitflag,output,lambda,grad,hessian]=fmincon('minCostimp6fit',x0,A,B,Aeq,Beq,LB,UB,'restrCost2imp6fit',options);
%x=x0;
%Cost structure
for k=1:n
Tx(k) =density(k)*[[1 x(k) x(k)^2 x(k)^3 x(k)^4 x(k)^5 x(k)^6] ...
     [1 x(k) x(k)^2 x(k)^3 x(k)^4 x(k)^5 x(k)^6]*(x(k+n)) ...
     [1 x(k) x(k)^2 x(k)^3 x(k)^4 x(k)^5 x(k)^6]*(x(k+n))^2 ...
     [1 x(k) x(k)^2 x(k)^3 x(k)^4 x(k)^5 x(k)^6]*(x(k+n))^3 ...
     [1 x(k) x(k)^2 x(k)^3 x(k)^4 x(k)^5 x(k)^6]*(x(k+n))^4 ...
     [1 x(k) x(k)^2 x(k)^3 x(k)^4 x(k)^5 x(k)^6]*(x(k+n))^5 ...
     [1 x(k) x(k)^2 x(k)^3 x(k)^4 x(k)^5 x(k)^6]*(x(k+n))^6]*theta; 
NL(k)=density(k)*(bw)/10;
Irayo(k)=31*((0.6*NL(k))/Tx(k)-1)^(1/2.6);
Rb(k)=sqrt(rho(k)*Eog*x(k)^2/(rho(k)*Eog-2*pi*Irayo(k)*x(k)^2));
Ig(k)=rho(k)*Eog/(2*pi*Rb(k)^2);
CG(k)=NT(k)*b(1)*(rho(k)/Rb(k))^(b(2));
CI(k)=NT(k)*NPH(k)*kI*(x(k+n)-w60);
% dTdR(k)=(density(k))*((a(2)+a(7)*(x(k+n))+a(12)*(x(k+n))^2+a(17)*(x(k+n))^3+a(22)*(x(k+n))^4)+...
%                      2*(a(3)+a(8)*(x(k+n))+a(13)*(x(k+n))^2+a(18)*(x(k+n))^3+a(23)*(x(k+n))^4)*x(k)+...
%                      3*(a(4)+a(9)*(x(k+n))+a(14)*(x(k+n))^2+a(19)*(x(k+n))^3+a(24)*(x(k+n))^4)*x(k)^2+...
%                      4*(a(5)+a(10)*(x(k+n))+a(15)*(x(k+n))^2+a(20)*(x(k+n))^3+a(25)*(x(k+n))^4)*x(k)^3);
% A1(k)=dTdR(k)*31*.6*NL(k)/Tx(k)^2;
% A2(k)=31*(1/2.6)*(0.6*NL(k)/Tx(k)-1)^(inv(2.6)-1)*A1(k);
% A3(k)=(1-2*pi*Irayo(k)*x(k)^2/(rho(k)*Eog))^-.5-.5*x(k)*(2*pi/(rho(k)*Eog))*(A2(k)*x(k)^2+2*x(k)*Irayo(k))*(1-2*pi*Irayo(k)*x(k)^2/(rho(k)*Eog))^-1.5;
% dCGdR(k)=-NT(k)*b(2)*b(1)*((rho(k)/Rb(k))^(b(2))-1)*rho(k)*inv(Rb(k)^2)*A3(k);
% dCIdW(k)=kI*NT(k)*NPH(k);
% dTdW(k)=density(k)*((a(6)+a(7)*x(k)+a(8)*x(k)^2+a(9)*x(k)^3+a(10)*x(k)^4)+...
%      2*(a(11)+a(12)*x(k)+a(13)*x(k)^2+a(14)*x(k)^3+a(15)*x(k)^4)*(x(k+n))^1+...
%      3*(a(16)+a(17)*x(k)+a(18)*x(k)^2+a(19)*x(k)^3+a(20)*x(k)^4)*(x(k+n))^2+...
%      4*(a(21)+a(22)*x(k)+a(23)*x(k)^2+a(24)*x(k)^3+a(25)*x(k)^4)*(x(k+n))^3);
% A4(k)=dTdW(k)*31*.6*NL(k)/Tx(k)^2;
% A5(k)=31*(1/2.6)*(0.6*NL(k)/Tx(k)-1)^(inv(2.6)-1)*A4(k);
% A6(k)=-.5*x(k)*((1-2*pi*Irayo(k)*x(k)^2/(rho(k)*Eog))^-1.5)*A5(k)*(2*pi*x(k)^2)/(rho(k)*Eog);
% dCGdW(k)=-NT(k)*b(2)*b(1)*((rho(k)/Rb(k))^(b(2))-1)*rho(k)*inv(Rb(k)^2)*A6(k);
% dCdT(k)=dCGdR(k)*(dTdR(k)).^-1;
count(k)=k;
% fl(u)=exitflag;
% Incremental(u,k)=dCdT;
% end
% Tasa(u)=sum(Tx.*L/Lt);
% CGt(u)=sum(CG);
end
%plot(CGt,Tasa,'r-')
%plot(Tasa,Incremental,'k-')
% xf=x;
% 
plot(count,CG+CI,'b-',count,CG,'g-',count,CI,'c--')
hold on
yyaxis right
plot(count,Tx,'r-',count,Tmax,'k-')

% for r=1:60;
% y0=200;
% y=fsolve('tasa',y0)
% wxx(r)=y;
% count(r)=r;
% Tr(r)=density(1)*(...
%     (a(1)+a(2)*r+a(3)*r^2+a(4)*r^3+a(5)*r^4)+...
%     (a(6)+a(7)*r+a(8)*r^2+a(9)*r^3+a(10)*r^4)*wxx(r)+...
%     (a(11)+a(12)*r+a(13)*r^2+a(14)*r^3+a(15)*r^4)*wxx(r)^2+...
%     (a(16)+a(17)*r+a(18)*r^2+a(19)*r^3+a(20)*r^4)*wxx(r)^3+...
%     (a(21)+a(22)*r+a(23)*r^2+a(24)*r^3+a(25)*r^4)*wxx(r)^4);
% fo(r)=0*NT(k)*(b(1)+b(2)*rho(k)+b(3)*rho(k)^2)*r^(b(4)) + 1*NT(k)*NPH(k)*(kI*(wxx(r)-w60));
% end
% figure
% plot(fo);
% figure
% plot(count,Tr,'b-')
% figure
% plot(count,wxx,'b-')
% % x

%%Lagrange
xf=x;
mu=lambda.ineqnonlin';
lmd=lambda.eqnonlin;
%x0=horzcat(xf,mu,lmd(1));
% options2=optimset('MaxIter',30000000,'MaxFunEvals',10000000);
% x=fsolve('incrementIMP',x0,options2);
%x=x0;
% for k=1:n% Derivatives
%     aa=(a(1)+a(6)*(x(k+n))+a(11)*(x(k+n))^2+a(16)*(x(k+n))^3+a(21)*(x(k+n))^4);
% ab=(a(2)+a(7)*(x(k+n))+a(12)*(x(k+n))^2+a(17)*(x(k+n))^3+a(22)*(x(k+n))^4);
% ac=(a(3)+a(8)*(x(k+n))+a(13)*(x(k+n))^2+a(18)*(x(k+n))^3+a(23)*(x(k+n))^4);
% ad=(a(4)+a(9)*(x(k+n))+a(14)*(x(k+n))^2+a(19)*(x(k+n))^3+a(24)*(x(k+n))^4);
% ae=(a(5)+a(10)*(x(k+n))+a(15)*(x(k+n))^2+a(20)*(x(k+n))^3+a(25)*(x(k+n))^4);
% Tx(k)=(density(k))*((a(1)+a(6)*(x(k+n))+a(11)*(x(k+n))^2+a(16)*(x(k+n))^3+a(21)*(x(k+n))^4)+...
%                     (a(2)+a(7)*(x(k+n))+a(12)*(x(k+n))^2+a(17)*(x(k+n))^3+a(22)*(x(k+n))^4)*x(k)+...
%                     (a(3)+a(8)*(x(k+n))+a(13)*(x(k+n))^2+a(18)*(x(k+n))^3+a(23)*(x(k+n))^4)*x(k)^2+...
%                     (a(4)+a(9)*(x(k+n))+a(14)*(x(k+n))^2+a(19)*(x(k+n))^3+a(24)*(x(k+n))^4)*x(k)^3+...
%                     (a(5)+a(10)*(x(k+n))+a(15)*(x(k+n))^2+a(20)*(x(k+n))^3+a(25)*(x(k+n))^4)*x(k)^4);     
% NL(k)=density(k)*(bw)/10;
% Irayo(k)=31*((0.6*NL(k))/Tx(k)-1)^(1/2.6);
% Rb(k)=sqrt((rho(k)*Eog*x(k)^2)/(rho(k)*Eog-2*pi*Irayo(k)*x(k)^2));   
% CG(k)=NT(k)*b(1)*(rho(k)/Rb(k))^(b(2));
% CI(k)=NT(k)*NPH(k)*kI*(x(k+n)-w60);
% dTdR(k)=(density(k))*((a(2)+a(7)*(x(k+n))+a(12)*(x(k+n))^2+a(17)*(x(k+n))^3+a(22)*(x(k+n))^4)+...
%                      2*(a(3)+a(8)*(x(k+n))+a(13)*(x(k+n))^2+a(18)*(x(k+n))^3+a(23)*(x(k+n))^4)*x(k)+...
%                      3*(a(4)+a(9)*(x(k+n))+a(14)*(x(k+n))^2+a(19)*(x(k+n))^3+a(24)*(x(k+n))^4)*x(k)^2+...
%                      4*(a(5)+a(10)*(x(k+n))+a(15)*(x(k+n))^2+a(20)*(x(k+n))^3+a(25)*(x(k+n))^4)*x(k)^3);
% % A1(k)=dTdR(k)*.6*NL(k)/Tx(k)^2;
% % A2(k)=31*(1/2.6)*(0.6*NL(k)/Tx(k)-1)^(inv(2.6)-1)*A1(k);
% % Ax(k)=A2(k)*x(k)^2+2*Irayo(k)*x(k);
% %  A3(k)=(Eog*rho(k)*x(k)*(pi*A2(k)*x(k)^3 + Eog*rho(k)))/((Eog*rho(k) -...
% %      2*pi*Irayo(k)*x(k)^2)^2*((Eog*rho(k)*x(k)^2)/(Eog*rho(k) - 2*pi*Irayo(k)*x(k)^2))^(1/2));
% % A3(k)=(1-(2*pi*Irayo(k)*x(k)^2)/(rho(k)*Eog))^(-1/2)+...
% %     0.5*x(k)*((1-(2*pi*Irayo(k)*x(k)^2)/(rho(k)*Eog))^-(3/2))*((2*pi)/(rho(k)*Eog)*Ax(k));
% % dCGdRx(k)=-NT(k)*b(2)*b(1)*((rho(k)/Rb(k))^(b(2)-1))*rho(k)*inv(Rb(k)^2)*A3(k); 
% 
% dCGdR(k)=-(Eog*NT(k)*b(1)*b(2)*rho(k)*x(k)*(rho(k)/((Eog*rho(k)*x(k)^2)/...
%     (Eog*rho(k) - 62*pi*x(k)^2*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x(k) +...
%     5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/...
%     (5*density(k)*(aa + ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 + ab*x(k))))^(5/13)))^(1/2))^b(2)*(13*...
%     Eog*aa^2*density(k)*rho(k)*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x(k) +...
%     5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*...
%     (aa + ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 + ab*x(k))))^(8/13) - 186*pi*NL(k)*ac*x(k)^4 -...
%     279*pi*NL(k)*ad*x(k)^5 - 372*pi*NL(k)*ae*x(k)^6 - 93*pi*NL(k)*ab*x(k)^3 +...
%     13*Eog*ab^2*density(k)*rho(k)*x(k)^2*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x(k) +...
%     5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*(aa +...
%     ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 + ab*x(k))))^(8/13) + 13*Eog*ac^2*density(k)*rho(k)*x(k)^4*(-(5*aa*density(k) -...
%     3*NL(k) + 5*ab*density(k)*x(k) + 5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*(aa +...
%     ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 + ab*x(k))))^(8/13) +...
%     13*Eog*ad^2*density(k)*rho(k)*x(k)^6*(-(5*aa*density(k) - 3*NL(k) +...
%     5*ab*density(k)*x(k) + 5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 +...
%     5*ae*density(k)*x(k)^4)/(5*density(k)*(aa + ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 +...
%     ab*x(k))))^(8/13) + 13*Eog*ae^2*density(k)*rho(k)*x(k)^8*(-(5*aa*density(k) - 3*NL(k) +...
%     5*ab*density(k)*x(k) + 5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 +...
%     5*ae*density(k)*x(k)^4)/(5*density(k)*(aa + ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 + ab*x(k))))^(8/13) +...
%     26*Eog*aa*ac*density(k)*rho(k)*x(k)^2*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x(k) +...
%     5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*(aa +...
%     ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 + ab*x(k))))^(8/13) + 26*Eog*aa*ad*density(k)*rho(k)*x(k)^3*(-(5*aa*density(k) -...
%     3*NL(k) + 5*ab*density(k)*x(k) + 5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*(aa +...
%     ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 + ab*x(k))))^(8/13) + 26*Eog*ab*ac*density(k)*rho(k)*x(k)^3*(-(5*aa*density(k) - 3*NL(k) +...
%     5*ab*density(k)*x(k) + 5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*(aa +...
%     ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 + ab*x(k))))^(8/13) + 26*Eog*aa*ae*density(k)*rho(k)*x(k)^4*(-(5*aa*density(k) - 3*NL(k) +...
%     5*ab*density(k)*x(k) + 5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*(aa + ac*x(k)^2 + ad*x(k)^3 +...
%     ae*x(k)^4 + ab*x(k))))^(8/13) + 26*Eog*ab*ad*density(k)*rho(k)*x(k)^4*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x(k) +...
%     5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*(aa + ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 +...
%     ab*x(k))))^(8/13) + 26*Eog*ab*ae*density(k)*rho(k)*x(k)^5*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x(k) +...
%     5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*(aa + ac*x(k)^2 + ad*x(k)^3 +...
%     ae*x(k)^4 + ab*x(k))))^(8/13) + 26*Eog*ac*ad*density(k)*rho(k)*x(k)^5*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x(k) +...
%     5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*(aa + ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 +...
%     ab*x(k))))^(8/13) + 26*Eog*ac*ae*density(k)*rho(k)*x(k)^6*(-(5*aa*density(k) - 3*NL(k) +...
%     5*ab*density(k)*x(k) + 5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*(aa + ac*x(k)^2 +...
%     ad*x(k)^3 + ae*x(k)^4 + ab*x(k))))^(8/13) + 26*Eog*ad*ae*density(k)*rho(k)*x(k)^7*(-(5*aa*density(k) - 3*NL(k) +...
%     5*ab*density(k)*x(k) + 5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*(aa +...
%     ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 + ab*x(k))))^(8/13) + 26*Eog*aa*ab*density(k)*rho(k)*x(k)*(-(5*aa*density(k) - 3*NL(k) +...
%     5*ab*density(k)*x(k) + 5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*(aa + ac*x(k)^2 +...
%     ad*x(k)^3 + ae*x(k)^4 + ab*x(k))))^(8/13)))/(13*density(k)*(Eog*rho(k) - 62*pi*x(k)^2*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x(k) +...
%     5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 + 5*ae*density(k)*x(k)^4)/(5*density(k)*(aa + ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 +...
%     ab*x(k))))^(5/13))^2*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x(k) + 5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 +...
%     5*ae*density(k)*x(k)^4)/(5*density(k)*(aa + ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 + ab*x(k))))^(8/13)*((Eog*rho(k)*x(k)^2)/(Eog*rho(k) -...
%     62*pi*x(k)^2*(-(5*aa*density(k) - 3*NL(k) + 5*ab*density(k)*x(k) + 5*ac*density(k)*x(k)^2 + 5*ad*density(k)*x(k)^3 +...
%     5*ae*density(k)*x(k)^4)/(5*density(k)*(aa + ac*x(k)^2 + ad*x(k)^3 + ae*x(k)^4 + ab*x(k))))^(5/13)))*(aa + ac*x(k)^2 +...
%     ad*x(k)^3 + ae*x(k)^4 + ab*x(k))^2);
% 
% dCIdW(k)=kI*NT(k)*NPH(k);
%  dTdW(k)=density(k)*((a(6)+a(7)*x(k)+a(8)*x(k)^2+a(9)*x(k)^3+a(10)*x(k)^4)+...
%       2*(a(11)+a(12)*x(k)+a(13)*x(k)^2+a(14)*x(k)^3+a(15)*x(k)^4)*(x(k+n))^1+...
%       3*(a(16)+a(17)*x(k)+a(18)*x(k)^2+a(19)*x(k)^3+a(20)*x(k)^4)*(x(k+n))^2+...
%       4*(a(21)+a(22)*x(k)+a(23)*x(k)^2+a(24)*x(k)^3+a(25)*x(k)^4)*(x(k+n))^3);
% A4(k)=dTdW(k)*.6*NL(k)/Tx(k)^2;
% A5(k)=31*(1/2.6)*(0.6*NL(k)/Tx(k)-1)^(inv(2.6)-1)*A4(k);
% A6(k)=-.5*x(k)*((1-2*pi*Irayo(k)*x(k)^2/(rho(k)*Eog))^-1.5)*A5(k)*(2*pi*x(k)^2)/(rho(k)*Eog);
% dCGdW(k)=-NT(k)*b(2)*b(1)*((rho(k)/Rb(k))^(b(2)-1))*rho(k)*inv(Rb(k)^2)*A6(k);
% end
% for k=1:n
% F(k)=dCGdR(k)+x(6*n+1)*(L(k)/Lt)*dTdR(k)+x(k+4*n)-x(k+5*n);
% F(k+n)=dCGdW(k)+dCIdW(k)+x(6*n+1)*(L(k)/Lt)*dTdW(k)-x(k+2*n)+x(k+3*n);
% F(k+2*n)=x(k+2*n)*(w60-(x(k+n)));
% F(k+3*n)=x(k+3*n)*((x(k+n))-wmax);
% F(k+4*n)=x(k+4*n)*(x(k)-R0);
% F(k+5*n)=x(k+5*n)*(0-x(k));
% end
T=0;
for k=1:n
T=T+Tx(k)*L(k)/Lt;
end
disp('Total cost:')
Ctot=sum(CG)+sum(CI)
disp('Total grounding cost')
sum(CG)
disp('Total insulation cost')
sum(CI)
disp('Impulse resistance (ohms) + Insulation length (cm) per section')
x
disp('BFOR (outages/100k/yr per section)')
Tx
