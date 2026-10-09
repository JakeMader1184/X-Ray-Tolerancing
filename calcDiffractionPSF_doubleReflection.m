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
% r0 = 1500e-3;
r0 = 250e-3;


stressMag = -0.1*0;
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



Ephoton = 1e3; % Photon energy [eV]
wDetector = 200e-6; % Total width of the detector [m]
zDetector = 0; % Z-coordinate of detector [m]

alpha = 1/4*atan(r0/z0);
deltaR1 = Lp*sin(alpha);
lambda = 4.135667662e-15*3e8/Ephoton; % Photon wavelength = hc/E [m]
Npoints1 = ceil(Lp/(lambda/4/pi/alpha/sin(alpha)));
Npoints2 = 2*ceil(Lp/(lambda*z0/4/pi/sin(alpha)/wDetector*2));
NpointsD = 2*ceil(wDetector/(lambda*z0/2/pi/sin(alpha)/Lp));

xDetector = linspace(-wDetector/2,wDetector/2,NpointsD);
% div=5;
% tidxVect = 2:div:ceil(size(surfP.t,2)/2);
tidxVect = 41;
PSF = zeros(length(xDetector),length(tidxVect));

for tidx = 1:length(tidxVect)
    z1 = linspace(min(surfP.z(:))+1e-7,max(surfP.z(:)),Npoints1);
    t1 = mean(surfP.t(:,tidxVect(tidx)))*ones(size(z1));
    x1 = surfP.r(t1,z1);
    
    z2 = linspace(min(surfS.z(:)),max(surfS.z(:))-1e-7,Npoints2);
    t2 = mean(surfS.t(:,tidxVect(tidx)))*ones(size(z2));
    x2 = surfS.r(t2,z2);
    
    
    
    tic

    E2 = zeros(size(z2));
    for idx = 1:length(z2)
        d12 = sqrt((x2(idx)-x1).^2+(z2(idx)-z1).^2);
        integrand = sqrt(x1./d12).*exp(-2*pi*1i/lambda*(d12-z1));
        E2(idx) = deltaR1/Lp/sqrt(lambda*x2(idx))*trapz(z1,integrand);
    end
    toc
    for idx = 1:length(xDetector)
        x = xDetector(idx);
        
        d20 = sqrt((x2-x).^2+(z2-zDetector).^2);
        
        integrand = E2.*sqrt(x2./d20).*exp(-2*pi*1i/lambda*(d20));
        PSF(idx,tidx) = deltaR1/Lp^2/lambda/r0*...
            abs(trapz(z2,integrand))^2/(180/pi*3600/z0); % [1/arcsec]
    end
    
end

%%
% Calculate HPD
% Find center of PSF
rVect = xDetector(xDetector>=0)/z0*180/pi*3600;
% [rVect, sortidx] = sort(abs(xDetector)/z0*180/pi*3600);
PSFall = mean(PSF,2);

EEF = cumtrapz(2*rVect,PSFall(xDetector>=0));
HPD = interp1(EEF,2*rVect,0.5);
% EEF = cumtrapz(2*rVect,PSFall(sortidx));
% [~,minidx] = min(abs(EEF-0.5));
% HPD = 2*rVect(minidx);
disp(HPD)

figure(1)
plot(xDetector/z0*180/pi*3600, PSF,'.-')
set(gca,'xlim',[-wDetector/2/z0*180/pi*3600,wDetector/2/z0*180/pi*3600]);

figure(2)
clf
plot(2*rVect, EEF,'.')
hold on
plot([0,HPD],[0.5,0.5],'k-')
plot([HPD,HPD],[0,0.5],'k-')
hold off
text(HPD*1.5,0.5,['HPD = ' num2str(HPD,'%4.3f') ' arc-sec'])

figure(3)
plot(xDetector/z0*180/pi*3600, PSFall,'.-')
set(gca,'xlim',[-wDetector/2/z0*180/pi*3600,wDetector/2/z0*180/pi*3600]);


% save(['Deformation data\' defFileName 'diffractionResults.mat'])
save(['Deformation data\' defFileName '_0.1keV_AperturediffractionResults.mat'])
