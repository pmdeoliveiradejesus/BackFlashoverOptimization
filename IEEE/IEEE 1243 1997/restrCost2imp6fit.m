function [C,Ceq]=myfun(x)
global theta Tmax density n L Lt Tx T w60 wmax R0
T=0;
for k=1:n
Tx(k) =density(k)*[[1 x(k) x(k)^2 x(k)^3 x(k)^4 x(k)^5 x(k)^6] ...
     [1 x(k) x(k)^2 x(k)^3 x(k)^4 x(k)^5 x(k)^6]*(x(k+n)) ...
     [1 x(k) x(k)^2 x(k)^3 x(k)^4 x(k)^5 x(k)^6]*(x(k+n))^2 ...
     [1 x(k) x(k)^2 x(k)^3 x(k)^4 x(k)^5 x(k)^6]*(x(k+n))^3 ...
     [1 x(k) x(k)^2 x(k)^3 x(k)^4 x(k)^5 x(k)^6]*(x(k+n))^4 ...
     [1 x(k) x(k)^2 x(k)^3 x(k)^4 x(k)^5 x(k)^6]*(x(k+n))^5 ...
     [1 x(k) x(k)^2 x(k)^3 x(k)^4 x(k)^5 x(k)^6]*(x(k+n))^6]*theta; 
end
for k=1:n
T=T+Tx(k)*L(k)/Lt;
end

 for k=1:n
 C(k+0*n)=-x(k)+0;
 end
 for k=1:n
 C(k+1*n)=-(x(k+n))+w60;
 end
%  
% 
%   
    for k=1:n
 C(k+2*n)=(x(k+n))-wmax;
  end
  
%   for k=1:n
%  C(k+3*n)=x(k)-R0;
%  end

% for k=1:n
% Ceq(n+k)=-Tmax+Tx(k);
% end
Ceq(1)=-Tmax+T;
% for k=1:n
% Ceq(1+k)=x(k+n)-263;
% end
% C=[];





%  for k=1:n
%  Ceq(k+n)=(x(k+n))-(wfix);
%  end
%  C=[];
% x
%  C
%   Ceq
%  pause
