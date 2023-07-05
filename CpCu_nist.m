function CpCu = CpCu_nist(T,m)
%Input: T in K and mass in kg
%Unit: JK^-1
% Temperature range: 4-300 K
% Source: NIST, https://trc.nist.gov/cryogenics/materials/OFHC%20Copper/OFHC_Copper_rev1.htm
% Note: Even if the NIST fit is valida between 4-300 K, it is a better approximation 
% to extend it to 1-400 K, and possibly even 1-1000 K.
   density = 8960; % [kg/m^3] No value in NIST database, so most trustworthy value
   if T<1
       T=1;
   end
   if (T<300)
      dc_a = -1.91844;
	  dc_b = -0.15973;
	  dc_c = 8.61013;
	  dc_d =-18.996; 
	  dc_e = 21.9661;
	  dc_f =-12.7328;
	  dc_g =  3.54322;
	  dc_h = -0.3797;
      logT=log10(T);
	  p = dc_a + dc_b .* (logT).^1 + dc_c .* (logT).^2  + dc_d .* (logT).^3 + dc_e.* (logT).^4 + dc_f .* (logT).^5 + dc_g .* (logT).^6 + dc_h .* (logT).^7;
      p = 10^p;
   else
      p = 361.5+0.093.*T;
   end
   CpCu= m.*p;	  
end