function [PinhD, domain] = Pinhib_noR_NDR(Rmax, Einhz, Epurz, dt, Ek)
% Uses saturated, Rmax repair rate only when R can't be resolved. 
% Einhz must be oriented with time as rows!!!
% P = relative photosynthesis due to UV inhibitions (0-1)
% dt = particle time step in seconds
% 2/29/12 R repair rate excluded. See Bwf_rmax_flat.m
% 5/14/12 Repair linearly down scaled for Epurz < Ek

PinhD = ones(size(Einhz));
domain = nan*ones(size(Einhz)); 
RmaxD = Rmax*ones(size(Epurz)); 
%rpr(Epurz<0.07) = 0; % finds where Epur too low for repair <0.07 W/m-2
frac = Epurz/Ek;
frac(Epurz>Ek) = 1;
RmaxD = frac.*RmaxD; %repair decreases for 
%RmaxD = rpr*Rmax; %for binary repair.

for j = 1:size(Einhz,2) %each particle
 P=1;
 
  for i = 1:size(Einhz,1) %time 
   if Einhz(i,j) <= Rmax && P==1 % if damage less than repair, no inhibition
    %P=P*1;
    domain(i,j) = 0; %no inhibition
   elseif Einhz(i,j) > Rmax   % if damage exceed repair, induce inhibition.
       P=P*(1-dt*Einhz(i,j)) + dt*RmaxD(i,j); %inhibition
       domain(i,j) = 1; %inhibit
   elseif Einhz(i,j) <= Rmax && P<1
       P = P+dt*RmaxD(i,j); %repair unless too dark
       domain(i,j) = 2; %
   else
       keyboard; % 3 inequalities above should cover all cases. 
       %domain(i,j) = 3;
   end
   
   if P>1 %correct for 'over repair' when dt*Rmax > 1-Pinh
       P = 1;
        domain(i,j) = domain(i,j)+0.5; 
   end
   
   PinhD(i,j) = P;

  end %time

end %particle



% Troubleshooting:
%    if P<0 %prevent P from going negative
%        PinhD(i,j) = 0; faili = [faili i]; failj = [failj j];
%    else
%        PinhD(i,j) = P;
%    end