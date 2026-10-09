clear all
% clc
addpath('..\Ray trace functions')

% Define geometric parameters, in addition to those defined in the
% deformation file
zgap = 0e-3; % Gap between primary and secondary mirrors
t0 = 0*pi/180; % The mirror must not cross +-pi, because this is where there is a discontinuity in the atan2 function

% If this is true, then the ray trace will be copied for an integer number
% of mirrors to complete a full shell.
replicateFullShellFlag = true;

% Integrated stress applied to front surface [N/m]
% These are magnitudes (peak values) of Legendre-Legendre polynomials up
% to 2nd order, with the ordering as: azimuth=[0,0,0,1,1,1,2,2,2] and
% axial=[0,1,2,0,1,2,0,1,2]
% stressMag = zeros(1,9);
% stressMag=-0.2;


hs = 0.5e-3;
Lp=100e-3;
stressFront=0.2;
r0 = 1500e-3;
% r0 = 250e-3;

stressMag = -0.1;
stressMagP = stressMag;
stressMagS = stressMag;

defFileName = ['T' num2str(round(hs*1e6)) 'R' num2str(r0*1000,'%04d') 'L' num2str(round(Lp*1e3)) 'Nf' num2str(abs(stressFront)*10,'%02d')];

% Load the deformation files and make sure the geometric parameters match
C = load(['Deformation data\' defFileName '.mat']);

Lp = C.geoparams.Lp;
r0 = C.geoparams.r0;
Wt = C.geoparams.Wt;
z0 = C.geoparams.z0;
Npointsz = 2*C.modelparams.nDivZ+1;
Npointst = 2*C.modelparams.nDivTheta+1;

% Define the rigid-body translation and rotation errors at the center of
% the primary and secondary mirrors. The translation error is [dr,dt,dz].
% The rotation error is [yaw, pitch, roll], and they are right-hand
% rotations about the positive [x,y,z] axes for t0=0.
% The angles are usually small so the order of rotation (which is yaw,
% pitch, then roll) does not have a significant effect.
eP = [0,0,0]*1e-3; % Translation error of the primary mirror [m]
phiP = [0,0,0]/3600*pi/180; % Rotation error of primary mirror [rad]
phiS = 3*phiP;
eS = eP;

% Create parameter structures
geoparams = struct('r0',r0,'t0',t0,'z0',z0,'Lp',Lp,'Wt',Wt,'zgap',zgap,...
    'Npointsz',Npointsz,'Npointst',Npointst,...
    'eP',eP,'eS',eS,'phiP',phiP,'phiS',phiS);
deformationparams = struct('defFileName',defFileName,'stressP',stressMagP,'stressS',stressMagS);

[surfP, surfS,zfocus] = GenMirrorSurface(geoparams,deformationparams);




deltaR = Lp*sin(1/4*atan(r0/z0));
Ephoton = 10e3; % Photon energy [eV]
lambda = 4.135667662e-15*3e8/Ephoton; % Photon wavelength = hc/E [m]

xDetector = linspace(-50e-6,50e-6,3000);
zDetector = -10.11188; % r0=1500e-3
% zDetector = -10.003; % r0=250e-3
zMirror = linspace(min(surfP.z(:)),max(surfP.z(:)),3000);
tMirror = zeros(size(zMirror));
xMirror = surfP.r(tMirror,zMirror);

PSF = zeros(size(xDetector));
for idx = 1:length(xDetector)
    x = xDetector(idx);
    d20 = sqrt((xMirror-x).^2+(zMirror-zDetector).^2);
    
    integrand = sqrt(xMirror./d20).*exp(-2*pi*1i/lambda*abs(d20-zMirror));
    PSF(idx) = deltaR/Lp^2/lambda/r0*...
        abs(trapz(zMirror,integrand))^2/(180/pi*3600/z0); % [1/arcsec]
end

% Calculate HPD
% Find center of PSF
[~,maxidx] = max(PSF);
[rVect, sortidx] = sort(abs((xDetector-xDetector(maxidx))/z0*180/pi*3600));

EEF = cumtrapz(2*rVect,PSF(sortidx));
[~,minidx] = min(abs(EEF-0.5));
HPD = 2*rVect(minidx);

figure(1)
plot(xDetector/z0*180/pi*3600, PSF,'.-')

figure(2)
clf
plot(2*rVect, EEF,'.')
hold on
plot([0,HPD],[0.5,0.5],'k-')
plot([HPD,HPD],[0,0.5],'k-')
hold off
text(HPD*1.5,0.5,['HPD = ' num2str(HPD,'%4.3f') ' arc-sec'])



