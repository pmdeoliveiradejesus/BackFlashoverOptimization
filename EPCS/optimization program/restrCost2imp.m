function [C,Ceq]=myfun(x)
global a Tmax density n L Lt Tx T w60 wmax R0
T=0;
for k=1:n
Tx(k)=(density(k))*((a(1)+a(6)*(x(k+n))+a(11)*(x(k+n))^2+a(16)*(x(k+n))^3+a(21)*(x(k+n))^4)+...
                    (a(2)+a(7)*(x(k+n))+a(12)*(x(k+n))^2+a(17)*(x(k+n))^3+a(22)*(x(k+n))^4)*x(k)+...
                    (a(3)+a(8)*(x(k+n))+a(13)*(x(k+n))^2+a(18)*(x(k+n))^3+a(23)*(x(k+n))^4)*x(k)^2+...
                    (a(4)+a(9)*(x(k+n))+a(14)*(x(k+n))^2+a(19)*(x(k+n))^3+a(24)*(x(k+n))^4)*x(k)^3+...
                    (a(5)+a(10)*(x(k+n))+a(15)*(x(k+n))^2+a(20)*(x(k+n))^3+a(25)*(x(k+n))^4)*x(k)^4);
end
for k=1:n
T=T+Tx(k)*L(k)/Lt;
end
 for k=1:n
 C(k+0*n)=-(x(k+n))+w60;
 end
 

  
    for k=1:n
 C(k+1*n)=(x(k+n))-wmax;
  end
  
  for k=1:n
 C(k+2*n)=x(k)-R0;
 end
 for k=1:n
 C(k+3*n)=-x(k)+0;
 end
% for k=1:n
% Ceq(n+k)=-Tmax+Tx(k);
% end
Ceq(1)=-Tmax+T;
% for k=1:175
% Ceq(1+k)=x(k+n)-263;
% end





%  for k=1:n
%  Ceq(k+n)=(x(k+n))-(wfix);
%  end
%  C=[];
% x
% C
% Ceq
% pause
