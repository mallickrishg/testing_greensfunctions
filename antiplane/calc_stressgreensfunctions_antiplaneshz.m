function [Stress_12,Stress_13] = calc_stressgreensfunctions_antiplaneshz(G,X,Y,x3_bot,L_x2,L_x3)
% OUTPUTS
% Stress_12 - stress components in response to e12 eigen strain 
%             contains sigma12,sigma13 as columns
% Stress_13 - stress components in response to e13 eigen strain 
%             contains sigma12,sigma13 as columns

% Boxcar function
boxc=@(x) (x+0.5>=0)-(x-0.5>=0);

% Stress kernels for distributed deformation (source(i)-receiver(j))
s1312 = @(D,L,W,x2,x3) G/(2*pi)*( ...
    log((x2-L/2).^2+(x3-D-W).^2) - log((x2+L/2).^2+(x3-D-W).^2) ...
    -log((x2-L/2).^2+(x3-D).^2) + log((x2+L/2).^2+(x3-D).^2));

s1212 = @(D,L,W,x2,x3) G/pi*( ...
    atan((x3-D)./(x2+L/2))-atan((x3-D)./(x2-L/2)) ...
    +atan((x3-D-W)./(x2-L/2))-atan((x3-D-W)./(x2+L/2))) ...
    - 2*G*boxc(x2/L).*boxc((x3-(2*D+W)/2)/W);

s1213 = @(D,L,W,x2,x3) G/(2*pi)*( ...
    log((x2-L/2).^2+(x3-D-W).^2) - log((x2+L/2).^2+(x3-D-W).^2) ...
    -log((x2-L/2).^2+(x3-D).^2) + log((x2+L/2).^2+(x3-D).^2));

s1313=@(D,L,W,x2,x3) G/pi*( ...
    atan((x2+L/2)./(x3-D))  -atan((x2-L/2)./(x3-D)) ...
    -atan((x2+L/2)./(x3-D-W))+atan((x2-L/2)./(x3-D-W))) ...
    - 2*G*boxc(x2/L).*boxc((x3-(2*D+W)/2)/W);

% calculate stresses from a shear zone for eigen strain in 12-direction
s12 = s1212(x3_bot,L_x2,L_x3,X,Y);
s13 = s1213(x3_bot,L_x2,L_x3,X,Y);
Stress_12 = [s12,s13];


% calculate stresses from a shear zone for eigen strain in 13-direction
s12 = s1312(x3_bot,L_x2,L_x3,X,Y);
s13 = s1313(x3_bot,L_x2,L_x3,X,Y);
Stress_13 = [s12,s13];

end