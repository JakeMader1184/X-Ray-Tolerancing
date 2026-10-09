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

r0Vect = [250:250:1500]*1e-3;

% Calculate effective area for each shell
etaVect = [0.931,0.865,0.802,0.739,0.676,0.612]; % Reflectivity for each shell
rlower = [250e-3, r0Vect(2:end)-125e-3]; % Radius halfway between each shell and the next smallest shell 
rupper = [r0Vect(1:end-1)+125e-3, 1500e-3]; % Radius halfway between each shell and the next largest shell
Ashell = pi*(rupper.^2-rlower.^2).*etaVect; % Relative effective area [m^2]

stressMagVect=-0.1*ones(size(r0Vect));

HPD = zeros(length(r0Vect),1);
RMSD = zeros(length(r0Vect),1);
drRMS = zeros(length(r0Vect),2);
drdzRMS = zeros(length(r0Vect),2);
drdtRMS = zeros(length(r0Vect),2);
rayAngle = [];
weight = zeros(size(r0Vect));

for idxR = 1:length(r0Vect)
    
    r0 = r0Vect(idxR);
    stressMag = stressMagVect(idxR);
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
    
    rayAngle = [rayAngle, reshape(2*abs(surfP.drdzFE-flipud(surfS.drdzFE)),[],1)*180/pi*3600];
    weight(idxR) = Ashell(idxR)/Npointsz/Npointst;
%     HPD(idxR) = 2*median(reshape(2*abs(surfP.drdzFE-flipud(surfS.drdzFE)),[],1))*180/pi*3600;
    HPD(idxR) = 2*median(rayAngle(:,idxR));
    RMSD(idxR) = 2*sqrt(mean(reshape(2*abs(surfP.drdzFE-flipud(surfS.drdzFE)),[],1).^2))*180/pi*3600;
    
    
    
    radfitP = pinv([ones(size(surfP.t(:))), cos(surfP.t(:))])*surfP.drFE(:);
    radfitS = pinv([ones(size(surfS.t(:))), cos(surfS.t(:))])*surfS.drFE(:);
    drPAxial = surfP.drFE-reshape([ones(size(surfP.t(:))), cos(surfP.t(:))]*radfitP,size(surfP.drFE));
    drSAxial = surfS.drFE-reshape([ones(size(surfS.t(:))), cos(surfS.t(:))]*radfitS,size(surfS.drFE));
%     drPAxial = surfP.drFE-mean(surfP.drFE); 
%     drSAxial = surfS.drFE-mean(surfS.drFE);
    
    drRMS(idxR,:) = [sqrt(mean(drPAxial(:).^2)), sqrt(mean(drSAxial(:).^2))];
%     drRMS(idxR,:) = [sqrt(mean(surfP.drFE(:).^2)), sqrt(mean(surfS.drFE(:).^2))];
    drdzRMS(idxR,:) = [sqrt(mean(surfP.drdzFE(:).^2)), sqrt(mean(surfS.drdzFE(:).^2))];
    drdtRMS(idxR,:) = [sqrt(mean(surfP.drdtFE(:).^2)), sqrt(mean(surfS.drdtFE(:).^2))];
       
        
    if idxR==1
        figure(4);
        surf(surfP.t*180/pi,(surfP.z-10.05)*1000,surfP.drFE*1e9,'edgecolor','none');
        set(gca,'fontsize',12,'ticklabelinterpreter','latex')
        view([0,90])
        xlabel('Azimuth angle [$^\circ$]','interpreter','latex')
        ylabel('Z [mm]','interpreter','latex')
        cbar = colorbar;
        cbar.Label.Interpreter='latex';
        cbar.Label.String = 'Radial displacement [nm]';
        cbar.Label.FontSize = 12;
        cbar.FontSize = 12;
        cbar.TickLabelInterpreter = 'latex';
        
        figure(5);
        surf(surfP.t*180/pi,(surfP.z-10.05)*1000,surfP.drdzFE*180/pi*3600,'edgecolor','none');
        set(gca,'fontsize',12,'ticklabelinterpreter','latex')
        view([0,90])
        xlabel('Azimuth angle [$^\circ$]','interpreter','latex')
        ylabel('Z [mm]','interpreter','latex')
        cbar = colorbar;
        cbar.Label.Interpreter='latex';
        cbar.Label.String = 'Axial slope [arc-sec]';
        cbar.Label.FontSize = 12;
        cbar.FontSize = 12;
        cbar.TickLabelInterpreter = 'latex';
    end
end
%%
tic
% Calculate total HPD
NrTest = 2000;
weightList = repmat(reshape(repmat(weight,size(rayAngle,1),1),[],1),1,NrTest);
rayAngleList = repmat(reshape(rayAngle,[],1),1,NrTest);
rAngleTestVect = linspace(0,max(rayAngle(:)),NrTest);
rAngleTest = rAngleTestVect;
% rAngleTest = repmat(linspace(0,max(rayAngle(:)),NrTest),size(rayAngleList,1),1);

EEF = sum(weightList.*(rayAngleList<=rAngleTest))/sum(weightList(:,1));
[minval, minidx] = min(abs(EEF-0.5));
HPDtotal = 2*rAngleTest(minidx);
toc


figure(1)
clf
plot(r0Vect*1000,HPD,'o-k','markerfacecolor','w')
colorOrder = get(gca,'colororder');
hold on
plot([0, max(r0Vect*1000)], [HPDtotal, HPDtotal],'--k')
hold off
set(gca,'fontsize',12,'ticklabelinterpreter','latex',...
    'xgrid','on','ygrid','on','xlim',[0,max(r0Vect*1000)],...
    'xtick',[0,r0Vect]*1000)
xlabel('Mirror radius [mm]','interpreter','latex')
ylabel('Half-power diameter [arc-sec]','interpreter','latex')
text(min(r0Vect*1000),HPDtotal*1.05,['Telescope HPD: ' num2str(HPDtotal,'%3.2f') ' arc-sec'],'interpreter','latex','verticalalignment','bottom')


figure(2)
clf
plot(r0Vect*1000,drRMS*1e9,'s-','markerfacecolor','w')
set(gca,'fontsize',12,'ticklabelinterpreter','latex',...
    'xgrid','on','ygrid','on','xlim',[0,max(r0Vect*1000)],...
    'xtick',[0,r0Vect]*1000)
xlabel('Mirror radius [mm]','interpreter','latex')
ylabel('RMS radial displacement [nm]','interpreter','latex')


figure(3)
clf
plot(r0Vect*1000,drdzRMS*180/pi*3600,'o-','markerfacecolor','w')
hold on
plot(r0Vect*1000,drdtRMS*180/pi*3600,'d-','markerfacecolor','w')
hold off
set(gca,'fontsize',12,'ticklabelinterpreter','latex',...
    'xgrid','on','ygrid','on','xlim',[0,max(r0Vect*1000)],...
    'xtick',[0,r0Vect]*1000)
xlabel('Mirror radius [mm]','interpreter','latex')
ylabel('RMS slope [arc-sec]','interpreter','latex')
leg = legend('Axial P', 'Axial S','Azimuthal P','Azimuthal S');
set(leg,'fontsize',12,'interpreter','latex','location','west')


