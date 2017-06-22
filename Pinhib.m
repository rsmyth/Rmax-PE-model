function [Pinh, domain] = Pinhib(Rmax, R, Einhz, dt)
% Uses saturated and unsaturated repair rates
% Einhz must be oriented with time as rows!!!
% P = relative photosynthesis due to UV inhibitions (0-1)
% dt = particle time step in seconds

Pinh = ones(size(Einhz));
domain = nan*ones(size(Einhz)); 
for j = 1:size(Einhz,2) %each particle
 P=1;
 
  for i = 1:size(Einhz,1) %time 
   if P > (1-Rmax/R)
% if repair is not saturated
    P=P*(1-dt*Einhz(i,j))+dt*R*(1-P);   %~eq. 7 from Fritz, etal. 2008
    domain(i,j) = 0;
   else
    P=P*(1-dt*Einhz(i,j)) + dt*Rmax;
   domain(i,j) = 1;
   end
   
  Pinh(i,j) = P;
%   P
%   dt*R*(1-P)
%   dt*Einhz(i,j)
  end %time

end %particle



% Troubleshooting:
%    if P<0 %prevent P from going negative
%        Pinh(i,j) = 0; faili = [faili i]; failj = [failj j];
%    else
%        Pinh(i,j) = P;
%    end