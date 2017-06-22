function [PUR Epar] = calcPUR(PA, Epi, PAR, z, alb, Kpar, dpar)

% Attenuates PAR with depth and calculates PUR correction as described in Neale et al. 2014
% Biogeosciences

%to check dimension shifts: 
%figure; plot(repmat(1-alb, 1,length(dpar)).*PAR); title('problem if * not on lines')


Qpi = Epi.*(dpar/119670);  % convert from mW/m2/nm to umol photon/m2/s/nm (=E^Q_pi)
deltad =gradient(dpar);
api = nansum(PA.*Qpi.*deltad)./nansum(Qpi.*deltad);  %eq. 5

%replicate into matrices of time x depth
alb = repmat(alb, [1 length(dpar)]);
dpar = repmat(dpar, [size(PAR,1) 1]);

%E just below surface
QPAR = (1-alb).*PAR.*(dpar/119670); % convert mW/m2/nm to umol photons/m2/s/nm
PAR = (1-alb).*PAR;

% replicate into matrices of wavelength x time x depth for computations below
deltad = repmat(deltad',[1 size(PAR,1) length(z)]);
PA = repmat(PA', [1 size(PAR,1) length(z)]); %=aplamda
Kpar = repmat(Kpar', [1 size(PAR,1) length(z)]);

PAR = repmat(PAR', [1 1 length(z)]);
QPAR = repmat(QPAR', [1 1 length(z)]);
z1 = repmat(z, [size(PAR,2) 1 size(PAR,1)]); z = shiftdim(z1, 2); 
 

Qz = QPAR.*exp(-Kpar.*z);    %umol photons/m2/s
Eparz = PAR.*exp(-Kpar.*z);   %mW/m2/nm

%hold on; plot(squeeze(Eparz(:,:,1))', '*'); 

ais = sum(PA.*Qz.*deltad)./sum(Qz.*deltad);  %eq. 6

%
PUR = squeeze(ais)/api;    %eq. 7. dimensionless fraction of utilizable PAR by time and depth              

Epar = squeeze(sum(Eparz/1000.*deltad)); %mW/m2 to W/m2 and sum over wavelengths for broadband PAR. 

 



