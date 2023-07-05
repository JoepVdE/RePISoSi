% ***** Self-inductance of a straight long rectangular bar ***** 
% 
% Reference: 
% [1] Piatek et al., "Self inductance of long conductor of rectangular 
%     cross section," Electrical Review, ISSN 0033-2097, R. 88 NR 8/2012. 
% 
% Inputs:
%   a = width (m) 
%   b = height (m) 
%   L = length (m) 
%
% *****
function [LSelf] = calcSelfInductanceRectangularBar(a, b, L)

    % Eq. (24) in [1]
    % LSelf = 2e-7*L*(log(2*L/(a + b)) + 13/12 ... 
    %     - 2/3*(b/a*atan(a/b) + a/b*atan(b/a)) ... 
    %     + 0.5*log(1 + a/b*2/(1 + (a/b)^2)) ... 
    %     + 1/12*((a/b)^2*log(1 + (b/a)^2) + (b/a)^2*log(1 + (a/b)^2)));

    
    % Faster version added by AnV on July 1, 2022
    c = a/b; 
    d = b/a; 
    % Eq. (24) in [1]
    LSelf = 2e-7*L*(log(2*L/(a + b)) + 13/12 ... 
        - 2/3*(d*atan(c) + c*atan(d)) ... 
        + 0.5*log(1 + c*2/(1 + (c)^2)) ... 
        + 1/12*((c)^2*log(1 + (d)^2) + (d)^2*log(1 + (c)^2)));
    
end
