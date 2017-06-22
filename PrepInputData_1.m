% Create input.mat for LES-BWF-PE model runs

clear

units =[];
 
% Select and load input data
pname = pwd;
cd([pname,'/Input/Irradiance'])

[Iname, Ipath, dex] = uigetfile('*.mat', 'Select a STAR model irradiance (mW/m2/nm) file');
units.Io = 'mW/m2/nm';
I = load([Ipath,'/',Iname]);  %(mW/m2/nm)    NOTE: Irradiance.mat time step uneven

cd([pname,'/LES'])

[Tname, Tpath, dex] = uigetfile('*.mat', 'Select a particle trajectory file');
p = load([Tpath,'/',Tname])

model = input('Enter 1 for Rmax with R or 2 for Rmax without R  ');
if model == 1
cd([pname,'/Input/BWF/RmaxR'])
elseif model == 2 
    cd([pname,'/Input/BWF/RmaxOnly'])
end
      
[Bname, Bpath, dex] = uigetfile('*.mat', 'Select a BWF file');
bwf = load([Bpath,'/',Bname]); 

if isfield(bwf, 'atn') %Kd appended to some BWFs
    atn = bwf.atn;
    
else
cd([pname,'/Input/Irradiance'])
[Aname, Apath, dex] = uigetfile('*.mat', 'Select an attenuation file');
atn = load([Apath,'/',Aname]); % attenuation coefs by wavelength (1/m) 
end

%if isfield(bwf, 'PA')    
if ~exist('PA') %PA appended to some BWFs
cd([pname,'/Input/Irradiance'])  
[PAname, PApath, dex] = uigetfile('*.mat', 'Select a PA file for PUR calculation');
PA = load([PApath,'/',PAname]); %for PUR calc
end

cd([pname,'/Input'])

doy = input('Enter start time in GMT Day of Year (e.g., 330.5 for R14B)  ');

% extract and transform input data

pcx = find(bwf.f_value >=6,1);  % find number of PC
epsilon = bwf.epsilon(pcx,:); units.eps = 'm2/J';             
epsPar = bwf.eps_par(pcx)/1000/1.18; % m2/J. PAR inhibition weight, 1.18, corrects for PAR range in the Beast (400-650). Added 2/14/12. Prior runs incorrect.  
Pbs = bwf.Ps(pcx)/3600; units.Pbs = 'gC/gChl/s';  %gC/gChl/h --> gC/gChl/s      
Ek = bwf.Ek(pcx)*1.18; units.Ek = 'J/s/m2';  %J/s/m2=W/m2 (aka Es)  1.18 corrects for PAR range in the Beast (400-650)
Rmax = bwf.Rmax;   % maximum repair rate /s
if model == 1
R = bwf.R_Emodel(pcx);   % repair rate /s
end
  
zp = p.z0;
% Kz trajectories are -z
if max(max(zp)) < 0  
    zp = -zp;
end

% Allign irradiance and particle trajectory data 
len = p.time(end)/86400; %length of trajectory file in fraction of day
ix = find(I.time >= doy & I.time<= doy+len);
time = p.time;  % seconds  %/86400+I.time(ix(1)); to convert particle time (s) to GMT of irradiance data 
d = I.d_mod;    % wavelength (nm)
z = [0:0.1:ceil(max(max(zp)))]; % depth for light field
I.Io_mod(isnan(I.Io_mod)) = 0;
Io = interp1(I.time, I.Io_mod, p.time/86400+I.time(ix(1))); % incident spectra interpolated to particle timestep    
fig = input('Enter 1 to see plot of irradiance data interpolated to model time step  ');
if fig == 1
    figure; plot(time/86400+I.time(ix(1)), Io,'-'); % hold on; plot((I.time(ix)-I.time(ix(1)))*86400, I.Io_mod(ix,[5:10:end]),'x');
end 
%clear I p

% albedo
ax = input('Enter 1 to use Fresnels eq to calculate albedo or 2 to enter a fixed albedo   ');
if ax == 2
alb = input('Enter fixed albedo (e.g, 0.05)  ');
alb = repmat(alb, length(p.time),1);
else
lat = input('Enter lattitude (e.g., -77 for R14B)   ');
[alb RefAng] = albedo(floor(doy)+(p.time/86400),lat);
end
fig = input('Enter 1 to see plot of albedo  ');
if fig == 1
    figure; plot(time/86400+I.time(ix(1)), alb,'-'); % hold on; plot((I.time(ix)-I.time(ix(1)))*86400, I.Io_mod(ix,[5:10:end]),'x');
end 


% Index UV and PAR wavelengths
px = find(d>=400); ux = find(d<=400);  %wavelength index for PAR and UV, respectively. NOTE: 400nm included in PAR and UV
if length(bwf.wave)> length(ux) %checks only for shorter wavelengths in BWF file. 
    bx = find(bwf.wave>=d(1)); 
else 
    bx = ux;
end
%clear bwf

% Extrapolate attenuation to all wavelengths
Kd = interp1(atn.d, atn.Kd,d,'pchip');   % = pchip(atn.d, atn.Kd,d);
i = find(gradient(Kd) <= 0 , 1); Kd(1:i) = Kd(i); % correct for decrease in Kd with pchip extrapolation to shorter wavelengths. 
fig = input('Enter 1 to see plot of interpolated Kd ');
if fig == 1
figure; plot(atn.d, atn.Kd, 'ro'); hold on; plot(d, Kd, 'k', 'linewidth', 1); xlabel('wavelength (nm)'); ylabel('Kd (1/m)');
end
%clear atn

source.irradiance = [Ipath Iname];
source.trajectories = [Tpath Tname];
source.bwf = [Bpath Bname];
if exist('Aname')
source.Kd = [Apath Aname];
end
if exist('PAname')
source.PA = [PApath PAname];
end

if model == 1
    outname = ['RmaxR_',Bname(1:3),'_',Tname(6:end-4)];
elseif model == 2
    outname = ['Rmax0_',Bname(1:3),'_',Tname(6:end-4)];
end

sv = input(['Save as ', outname,'.mat?  1 = yes, 2 = no, 3 = enter new file name  ']);

while sv == 3
    outname = input('Enter new name in single quotes  '  ); 
    sv = 1;
end

if sv == 1
    save(outname, 'source', 'Ek', 'Io','Kd', 'PA', 'Pbs', 'Rmax', 'alb', 'bx', 'd', 'doy', 'epsPar', 'epsilon',...
        'px','ux', 'time', 'units', 'zp');
else
    disp('File not saved')
end

cd(pname)
