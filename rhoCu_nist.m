% -------------------------------------------------------------------------
%  rhoCu_nist.m  --  NIST cryogenic resistivity of OFHC copper, including
%                    a magneto-resistance correction (Kohler-rule fit).
%  Part of RePISoSi - https://github.com/JoepVdE/RePISoSi  -  License: MIT
%  Author : J.L. Van den Eijnden, 2026
%           Original implementation by M. Mentink (Oct. 2022),
%           extended by J.L. Van den Eijnden and A. Vaskuri (Oct. 2023).
%
%  Inputs:
%      T    : temperature (K)         (scalar or vector)
%      B    : magnetic flux density (T) (scalar or vector)
%      RRR  : residual resistivity ratio (-) (scalar or vector)
%      f_MR : optional scaling factor for the magneto-resistivity (default 1)
%
%  Output:
%      rhoCu : resistivity (Ohm*m), broadcast to the largest input shape.
% -------------------------------------------------------------------------
function rhoCu = rhoCu_nist(T, B, RRR, f_MR)
    if nargin < 4
        f_MR = 1;   % default: no scaling of magneto-resistivity
    end

% % % If one of T, B, and RRR is given as a vector and one of two of the others are not, make them vectors with as many elements
if length(T)>=length(B) && length(T)>=length(RRR)
    rhoCu=zeros(size(T));

    if size(B)==[1,1]
        B=ones(size(T))*B;
    end
    if size(RRR)==[1,1]
        RRR=ones(size(T))*RRR;
    end
elseif length(B)>=length(T) && length(B)>=length(RRR)
    rhoCu=zeros(size(B));

    if size(RRR)==[1,1]
        RRR=ones(size(B))*RRR;
    end
    if size(T)==[1,1]
        T=ones(size(B))*T;
    end
elseif length(RRR)>=length(T) && length(RRR)>=length(B)
    rhoCu=zeros(size(RRR));

    if size(T)==[1,1]
        T=ones(size(RRR))*T;
    end
    if size(B)==[1,1]
        B=ones(size(RRR))*B;
    end
end

% % % Magnetic field direction is not important
B=abs(B);

% % % At very low magnetic field, the NIST fit causes numerical instability
idxLowB=find(B<0.1);
idxHighB=find(B>=0.1);

rho0=1.553e-8./RRR;
rhoi=1.171e-17*(T.^4.49)./(1+4.48e-7*(T.^3.35).*exp(-(50./T).^6.428));
rhoiref=0.4531*rho0.*rhoi./(rho0+rhoi);
rhcu=rho0+rhoi+rhoiref;

% Elements with very low magnetic field (B<0.1 T)
rhoCu(idxLowB)=rhcu(idxLowB);

% Other elements (B>=0.1 T)
lgs=0.43429*log(1.553D-8*B(idxHighB)./rhcu(idxHighB));
polys=-2.662+lgs.*(0.3168+lgs.*(0.6229+lgs.*(-0.1839+lgs*0.01827)));
corrs=(10.^polys);
rhoCu(idxHighB)=(1.+corrs*f_MR).*rhcu(idxHighB);

end


