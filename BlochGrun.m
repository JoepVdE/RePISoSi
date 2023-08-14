function rho = BlochGrun(x,A,thetaR,z,n)
%A=coeffs(1); 
%thetaR=coeffs(2); 
%n=coeffs(3); 
%s=coeffs(4); 
%z=coeffs(5);
%n=5;
%z=1E-12;
TLowLimit=2;
rhovalsout=zeros(size(x,2),size(x,1)); 




for i=1:length(x) %% Added by Matthias, 28/7/23, to get rid of annoying warnings
	if(x(i)<TLowLimit) 
		x(i)=TLowLimit;
	end
end

for i=1:length(x) 
       rhovalsout(i)=A.*(x(i)./thetaR).^n.*integral(@(x)((x.^n.*exp(x))./((exp(x)-1).^2)),0,thetaR./x(i))+z; 

end 
rho=rhovalsout';

end
