function F = myfun(x)
global a b kI rho n w60 NPH NT Tx Irayo density bw Eog NL k
F=0;
for k=1:n
 
Tx(k)=(density(k))*((a(1)+a(6)*(x(k+n))+a(11)*(x(k+n))^2+a(16)*(x(k+n))^3+a(21)*(x(k+n))^4)+...
                    (a(2)+a(7)*(x(k+n))+a(12)*(x(k+n))^2+a(17)*(x(k+n))^3+a(22)*(x(k+n))^4)*x(k)+...
                    (a(3)+a(8)*(x(k+n))+a(13)*(x(k+n))^2+a(18)*(x(k+n))^3+a(23)*(x(k+n))^4)*x(k)^2+...
                    (a(4)+a(9)*(x(k+n))+a(14)*(x(k+n))^2+a(19)*(x(k+n))^3+a(24)*(x(k+n))^4)*x(k)^3+...
                    (a(5)+a(10)*(x(k+n))+a(15)*(x(k+n))^2+a(20)*(x(k+n))^3+a(25)*(x(k+n))^4)*x(k)^4); 
NL(k)=density(k)*(bw)/10;
Irayo(k)=31*((0.6*NL(k))/Tx(k)-1)^(1/2.6);
Rb(k)=sqrt(rho(k)*Eog*x(k)^2/(rho(k)*Eog-2*pi*Irayo(k)*x(k)^2));   
CG(k)=NT(k)*b(1)*(rho(k)/Rb(k))^(b(2));
CI(k)=NT(k)*NPH(k)*kI*(x(k+n)-w60);
F=F+CG(k)+CI(k);
end

        
