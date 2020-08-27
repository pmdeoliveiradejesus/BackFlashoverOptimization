function F = myfun(x)
global theta b kI rho n w60 NPH NT Tx Irayo density bw Eog NL k
F=0;
%(density(k))*
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
CG(k)=NT(k)*b(1)*(rho(k)/Rb(k))^(b(2));
CI(k)=NT(k)*NPH(k)*kI*(x(k+n)-w60);
F=F+CG(k)+CI(k);
end

        
