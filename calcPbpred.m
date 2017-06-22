function [Pbpred] = calcPbpred(UV, d, bwf, kd, z, hr, PAR, dpar, alb, PUR)

% Input
% UV = UV irradiance at surface (W/m2 nm)
% d = UV wavelengths (nm)
% bwf = biological weighting function [i.e., epsilon(delta)] (m2/W)
% kd = attenuation of UV at each wavelength
% z = depth (m)
% hr = time (hrs)
% PAR = pohotosynthetically available radiation (W/m2 nm)
% alb = albedo (fraction of 1)
% PUR = PUR correction factors (in quanta)


bwf = repmat(bwf,[1,length(hr), length(z)]);
UV = repmat(UV, [1,1,length(z)]); 
kd = repmat(kd,[1,length(hr), length(z)]);
zrep  = repmat(z',[1, length(d), length(hr)]); zrep = shiftdim(zrep,1);
% PA = repmat(PA, [1,length(hr)]);
% Beast = repmat(Beast,[1,length(hr)]);


UVz = UV.*alb.*bwf.*exp(-kd.*zrep);  %UV irradiance at depth
Euv = sum(UVz,1); Euv = squeeze(Euv);

dx = find(dpar >= 651, 1);  %upper wavelength limit of Beast
for j = 1:length(z)  
Epar(:,j) = sum(PAR(1:dx,:)).*5./1000.*alb.*exp(-kpar.*(z(j)));  % wavelengths to 651
%Epur(:,j) = sum(PA.*PAR.*exp(-kpar.*z(j)));                     % Lehmann, et al. 2004 eq. 6 
end


Einh = Euv+(Epar*epsPar);
perinh = Euv*100./Einh;  %percent inhibition from UV


% Copied from xls for now.
PUR = [0.834190	0.830240	0.818700	0.801860	0.780720	0.759590	0.736750	0.713920	0.692640	0.671370	0.650090 0.568519185	0.516035828	0.483301114];
PUR = repmat(PUR, [length(hr),1]);
Pb = 1.7948978;
Ek = 4.96;


X = 1./Einh; X(X>1) = 1; 
Pbpred = Pb.*(1-exp(-Epar.*PUR./Ek)).*X;


%% Basic plots
%figure; plot(squeeze(UVz(:,5,:)) , -z)
figure; plot(hr, Pbpred(:,1:4)); hold on; plot(hr, Pbpred(:,5:8), '--'); plot(hr, Pbpred(:,9:11), ':');
legend('0m','1m','2m','3m','4m','5m','6m','7m','8m','9m','10m')
title('Predicted P^B'); xlabel('Hr');

figure; plot(Pbpred([1 3 5 7 9 11],:)', z); set(gca, 'ydir', 'rev'); legend('8', '10', '12','14','16','18') 
%plot(Pbpred(1,:)', z, 'k'); hold on; plot(Pbpred(2:6,:)', z); plot(Pbpred(7:11,:)', z, '--'); set(gca, 'ydir', 'rev'); 

figure; plot(hr,Euv); figure; plot(hr, Epar*epsPar); figure; plot(hr, Einh)


% Iz0 = UV.*alb.*exp(-kd.*zrep);
% hold on
% plot(squeeze(Iz0(20,5,:)) , -z, 'g') 