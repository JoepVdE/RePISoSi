% -------------------------------------------------------------------------
%  CpCu_nist.m  --  Heat capacity of OFHC copper (NIST cryogenic fit).
%  Part of RePISoSi - https://github.com/JoepVdE/RePISoSi  -  License: MIT
%  Author : J.L. Van den Eijnden, 2026
%           Original implementation by M. Mentink (Oct. 2022),
%           extended by J.L. Van den Eijnden and A. Vaskuri (Oct. 2023).
%
%  Inputs:
%      T : temperature (K), scalar
%      m : mass (kg), scalar
%  Output:
%      CpCu : heat capacity (J/K)
%
%  Source: NIST cryogenic materials database, OFHC Copper:
%  https://trc.nist.gov/cryogenics/materials/OFHC%20Copper/OFHC_Copper_rev1.htm
%  The NIST fit is formally valid 4-300 K; extending up to ~400 K (or beyond)
%  is a reasonable approximation in practice.
% -------------------------------------------------------------------------
function CpCu = CpCu_nist(T, m)
   density = 8960;  %#ok<NASGU>  % kg/m^3, kept for reference
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