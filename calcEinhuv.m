function [Einhuv] = calcEinhuv(BWF, Euv, time, alb, Kd, z, units)
% input data oriented such that time is by row and wavelengths by column

%replicate into matrices of wavelength x time x depth
%Einhuv = nan*ones(length(BWF), length(time), length(z));
Kd =   repmat(Kd',[1 length(time) length(z)]);  % 1/m
alb = repmat(alb', [length(BWF) 1 length(z)]);
BWF = repmat(BWF', [1 length(time) length(z)]); % Rmax(1/Jm2), T(1/Wm-2nm-1)
Euv = repmat(Euv', [1 1 length(z)])/1000;   % convert from mW/m2 to W/m2
z1 = repmat(z, [size(Euv,2) 1 size(Euv,1)]); z = shiftdim(z1, 2); 

if strcmp(units.eps, '1/(mWm-2nm-1)');
    BWF = BWF*1000; %1/mWm-2 --> 1/Wm-2
end

 
Einhuv = (1-alb).*BWF.*Euv.*exp(-Kd.*z);  %units = 1/(s nm) for BWF in m2/J
 

Einhuv = squeeze(nansum(Einhuv));  %deltad=1nm --> units = 1/s for Rmax, dimensionless for T 