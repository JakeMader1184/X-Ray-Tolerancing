function [FEdata, FEinterp]= calcConicalDeformation(geoparams, matparams, modelparams, mirrorPosition)
% Generates and runs a batch file for a conical substrate with a 
% film on the convex side of the substrate. The model is run with no stress
% applied to the film, to obtain the finite element nodal coordinates.
%
% INPUTS:
% geoparams is a structure containing geometric parameters
% matparams is a structure containing material parameters
% meshparams is a structure containing mesh parameters
% Calculates the deformation of a mirror resulting from film stress applied
% to the surfaces.
% 
% INPUTS:
% geoparams - a structure containing the following geometric parameters:
%   Type:   WI: Wolter type I
%           WS: Wolter-Schwartzchild [NOT CURRENTLY IMPLEMENTED]
%   z0:     Axial distance [m] between the junction of the primary and 
%           secondary mirrors and the nominal focal plane.
%   r0:     Radius [m] of junction of primary and secondary mirrors.
%   t0:     Azimuth angle [rad] of mirror centers
%   Lp:     Length [m] of primary mirror.
%   Wt:     Azimuthal extent [rad] of both mirrors
%   Npointsz: Number of points in the axial (Z) direction
%   Npointst: Number of points in the azimuthal (T) direction
%   eP, eS: Vectors of mirror translation errors [ex,ey,ez] [m] 
%   phiP, phiS: Vectors of mirror rotation errors [phix,phiy,phiz].
% matparams - a structure containing the following material parameters:
%   E:      Elastic modulus of material [Pa]
%   nu:     Poisson's ratio of material
% modelparams - a structure containing the following model parameters
%   nDivTheta, nDivZ: Number of elements in the azimuthal and axial
%                     directions. The number of nodes will be (2*n+1).
%   hs:     Substrate thickness [m]
%   hf:     Film thickness on front (convex) surface [m]
%   hb:     Film thickness on back (concave) surface [m]
%   constraintPos: a matrix (size Nx2) containing the positions, [tC, zC]
%               of each constraint. tC [rad] is the azimuthal coordinate,
%               and zC [m] is the axial coordinate.
%   constraintDOF: a character cell array (size Nx3) determining which degreees of 
%                  freedom at the corresponding contraint point are 
%                  fixed. Some examples:
%                       'XYZTXYZR' (='ALL'): All 6 DOFs fixed.
%                       'XYT':  X and Y translations fixed
%                       'ZR':   Z rotation fixed.
%   stressFront: integrated stress [N/m] applied to the front (concave)
%                mirror surface.
%   stressBack:  integrated stress [N/m] applied to the back (convex)
%                mirror surface. [NOT CURRENTLY IN USE]
%   realignFlag: Flag to determine whether the mirror should be re-aligned
%                to remove rigid-body displacements. This should typically
%                be true if the constraints are kinematic.
% 
% OUTPUTS:
% FEdata -  a structure containing the finite element data, with the
%           following fields:
%   nodes:  node numbers
%   x,y,z:  nodal coordinates [m]
%   dx,dy,dz: displacements [m]
%   drdz, drdt: rotations about the local element coordinates [rad]
% 
addpath('..\Ray trace functions') 
E = matparams.E;
nu = matparams.nu;
hs = modelparams.hs;
hf = modelparams.hf;
hb=0;

% Ef = E;
% nuf = nu;
Ef = 528e9;
nuf = 0.26;

% % Elastic properties of silicon in the (100) wafer orientation
% Exx = 169e9;
% Eyy = 169e9;
% Ezz = 130e9;
% nuxy=0.064;
% nuxz = 0.28;
% nuyz = 0.36;
% Gxy = 50.9e9;
% Gyz = 79.6e9;
% Gxz = 79.6e9;

% Set to produce -1 N/m equibiaxial stress per +1 K temperature change.
% For an isotropic elastic material, ADINA requires a CTE>=0.
CTEf = (1-nuf)/Ef/hf; 


% Generate mirror surfaces at the corners only.
geoparams.eP = [0 0 0]; geoparams.eS = [0 0 0]; 
geoparams.phiP = [0 0 0]; geoparams.phiS = [0 0 0];
geoparams.Npointsz = 10; geoparams.Npointst = 10;
deformationparams = struct('defFileName',[],'stressMag',0);

if mirrorPosition=='P'
    [surfData, ~,~] = GenMirrorSurface(geoparams, deformationparams);
elseif mirrorPosition=='S'
    [~,surfData,~] = GenMirrorSurface(geoparams, deformationparams);
end


% Calculate the in and out radii for each mirror
RinP = surfData.r(0,max(surfData.z(:)))+hs/2;
RoutP = surfData.r(0,min(surfData.z(:)))+hs/2;
% L = max(surfData.z(:))-min(surfData.z(:));
L = geoparams.Lp;

nDivTheta = modelparams.nDivTheta;
nDivZ = modelparams.nDivZ;

% Calculate the constraint point radii (on the substrate middle surface)
Nconstraints = size(modelparams.constraintPos,1);
constraintPosT = modelparams.constraintPos(:,1);
constraintPosZ = modelparams.constraintPos(:,2);
constraintPosR = (RinP-RoutP)/L*(constraintPosZ+L/2)+RoutP;
constraintPosRf = constraintPosR - hs/2 + hb/2;
constraintPosRb = constraintPosR + hs/2 - hf/2;
constraintLabel = 14+(1:Nconstraints)';
constraintLabelf = max(constraintLabel)+(1:Nconstraints)';
constraintLabelb = max(constraintLabelf)+(1:Nconstraints)';

constraintPos = [constraintLabel, constraintPosR, constraintPosT*180/pi, constraintPosZ];
constraintPosf = [constraintLabelf, constraintPosRf, constraintPosT*180/pi, constraintPosZ];
constraintPosb = [constraintLabelb, constraintPosRb, constraintPosT*180/pi, constraintPosZ];

batchStart = {
'database new save=no prompt=no';
'feprogram adina';
% *
'control undo=-1 automrebuild=no sessionstorage=no erroraction=skip';
'filesession no';
'control updatethickness=no';
'master idof=000000 initialstress=deformation epsi-first=yes';
'system name=1 type=cylindrical mode=1 AX=0 AY=0 AZ=1 BX=1 BY=0 BZ=0';
% *
'printout input-data=1 displacements=yes velocities=no accelerations=no';
'print-steps';
'entries block first last increment';
'1 1 1 1';
'dataend';

% *
% * Define geometry
% * ALL DIMENSIONS IN METERS
% *
'coordinates point system=1';
'entries name r theta xl';
% * Axis points
['1	0	0	' num2str(L/2)];
['2	0	0	' num2str(-L/2)];
% * Substrate middle surface points
['3	' num2str(RinP,'%10.10f') ' ' num2str(-geoparams.Wt/2*180/pi,'%10.10f') ' ' num2str(L/2,'%10.10f')];
['4	' num2str(RinP,'%10.10f') ' ' num2str(geoparams.Wt/2*180/pi,'%10.10f') ' ' num2str(L/2,'%10.10f')];
['5	' num2str(RoutP,'%10.10f') ' ' num2str(geoparams.Wt/2*180/pi,'%10.10f') ' ' num2str(-L/2,'%10.10f')];
['6	' num2str(RoutP,'%10.10f') ' ' num2str(-geoparams.Wt/2*180/pi,'%10.10f') ' ' num2str(-L/2,'%10.10f')];
% * Front film middle surface points
['7	' num2str(RinP-hs/2+hb/2,'%10.10f') ' ' num2str(-geoparams.Wt/2*180/pi,'%10.10f') ' ' num2str(L/2,'%10.10f')];
['8	' num2str(RinP-hs/2+hb/2,'%10.10f') ' ' num2str(geoparams.Wt/2*180/pi,'%10.10f') ' ' num2str(L/2,'%10.10f')];
['9	' num2str(RoutP-hs/2+hb/2,'%10.10f') ' ' num2str(geoparams.Wt/2*180/pi,'%10.10f') ' ' num2str(-L/2,'%10.10f')];
['10	' num2str(RoutP-hs/2+hb/2,'%10.10f') ' ' num2str(-geoparams.Wt/2*180/pi,'%10.10f') ' ' num2str(-L/2,'%10.10f')];
% % * Back film middle surface points
% ['11	' num2str(RinP+hs/2-hf/2,'%10.10f') ' ' num2str(-geoparams.Wt/2*180/pi,'%10.10f') ' ' num2str(L/2,'%10.10f')];
% ['12	' num2str(RinP+hs/2-hf/2,'%10.10f') ' ' num2str(geoparams.Wt/2*180/pi,'%10.10f') ' ' num2str(L/2,'%10.10f')];
% ['13	' num2str(RoutP+hs/2-hf/2,'%10.10f') ' ' num2str(geoparams.Wt/2*180/pi,'%10.10f') ' ' num2str(-L/2,'%10.10f')];
% ['14	' num2str(RoutP+hs/2-hf/2,'%10.10f') ' ' num2str(-geoparams.Wt/2*180/pi,'%10.10f') ' ' num2str(-L/2,'%10.10f')];
'dataend';

% *
% * Substrate lines
'line arc name=1 mode=1 p1=3 p2=4 center=1';
'line straight name=2 4 5';
'line arc name=3 mode=1 p1=5 p2=6 center=2';
'line straight name=4 6 3';

% * Film lines
'line arc name=5 mode=1 p1=7 p2=8 center=1';
'line straight name=6 8 9';
'line arc name=7 mode=1 p1=9 p2=10 center=2';
'line straight name=8 10 7';

% 'line arc name=9 mode=1 p1=11 p2=12 center=1';
% 'line straight name=10 12 13';
% 'line arc name=11 mode=1 p1=13 p2=14 center=2';
% 'line straight name=12 14 11';

% * Substrate mid-surface
'surface patch name=1 1 2 3 4';
% * Film surfaces
'surface patch name=2 5 6 7 8';
% 'surface patch name=3 9 10 11 12';

% *
'sfthickness';
'entries name thick';
['1 ' num2str(hs-hf-hb,'%10.10f')];
['2 ' num2str(hf,'%10.10f')];
% ['3 ' num2str(hb,'%10.10f')];
'dataend';
% *
'printnodes surfaces';
'entries surface';
'1';
'dataend';

% * Define coordinate axes for material
'axes constant name=1 ax=0 ay=1 az=1 bx=0 by=0 bz=1';
'set-axes-material surfaces';
'entries surface axes adir bdir';
'1 1 1 2';
'dataend';

% *
% * Define materials
% *
% ['material orthotropic name=1 ea=' num2str(Exx,'%10.10f') ' eb=' num2str(Eyy,'%10.10f') ' ec=' num2str(Ezz,'%10.10f') ' nuab=' num2str(nuxy,'%10.10f') ' nuac=' num2str(nuxz,'%10.10f') ' nubc=' num2str(nuyz,'%10.10f') ' Gab=' num2str(Gxy,'%10.10f') ' Gac=' num2str(Gxz,'%10.10f') ' Gbc=' num2str(Gyz,'%10.10f')];
% *
% * Define materials
% *
['material elastic name=1 e=' num2str(E,'%10.10f') ' nu=' num2str(nu,'%10.10f') ' alpha=' num2str(0,'%10.10f')];
['material elastic name=2 e=' num2str(Ef,'%10.10f') ' nu=' num2str(nuf,'%10.10f') ' alpha=' num2str(CTEf,'%10.10f')];
% ['material elastic name=3 e=' num2str(E,'%10.10f') ' nu=' num2str(nu,'%10.10f') ' alpha=' num2str(CTEb,'%10.10f')];


% *
% * Define element groups
% *
'kinematics displacements=small';
['egroup shell name=1 material=1 thickness=' num2str(hs-hf-hb,'%10.10f')];
['egroup shell name=2 material=2 thickness=' num2str(hf,'%10.10f') ' initialstrain=nodal stressreference=local'];
% ['egroup shell name=3 material=3 thickness=' num2str(hb,'%10.10f') ' initialstrain=nodal stressreference=local'];

% *
% * Subdivide surfaces
% *
['subdivide surface name=1 mode=divisions ndiv1=' num2str(nDivTheta) ' ndiv2=' num2str(nDivZ)];
'1';
'2';
% '3';
'dataend';

% *
% * Mesh surfaces
% *
'gsurface name=1 nodes=9 group=1';
'gsurface name=2 nodes=9 group=2';
% 'gsurface name=3 nodes=9 group=3';

% *
% * Make rigid links between surfaces
% *
'rigidlink name=1 slavetype=surface slavename=2 mastertype=surface mastername=1 link-opt=closest';
% 'rigidlink name=2 slavetype=surface slavename=3 mastertype=surface mastername=1 link-opt=closest';

% *
% * Apply temperature load to film
% *     temperature should be opposite of the stress value, since CTE>=0.
['load temperature name=1 magnitude=' num2str(-modelparams.stressFront)];
'apply-load';
'entries name ltype lname stype sname';
'1 temperature 1 surface 2';
% '2 temperature 2 surface 3';
'dataend';
% * Create boundary conditions
'fixity XT';
'''x-translation''';
'dataend';

'fixity XYT';
'''x-translation''';
'''y-translation''';
'dataend';

'fixity XYZT';
'''x-translation''';
'''y-translation''';
'''z-translation''';
'dataend'

'fixity XYZTZR';
'''x-translation''';
'''y-translation''';
'''z-translation''';
'''z-rotation''';
'dataend'

% *
% * Define constraint points
% *
'coordinates point system=1';
'entries name r theta xl'};

% *
% * Save file
% *
batchEnd = {
'dataend';
'adina file=''ConicalSubstrate.dat''';
'database save permfile=''ConicalSubstrate.idb'''};

fid = fopen('ConicalSubstrate.in','W');
fprintf(fid,'%s\r\n',batchStart{:});
fprintf(fid,'%d %f %f %f\r\n',[constraintPos; constraintPosf; constraintPosb]');
fprintf(fid,'%s\r\n','dataend', 'node-snap group=1');
fprintf(fid,'%d\r\n',constraintLabel');
fprintf(fid,'%s\r\n','dataend', 'node-snap group=2');
fprintf(fid,'%d\r\n',constraintLabelf');
% fprintf(fid,'%s\r\n','dataend', 'node-snap group=3');
% fprintf(fid,'%d\r\n',constraintLabelb');
% For each constraint point, apply the appropriate fixity
fprintf(fid,'%s\r\n','dataend','fixboundary points XT');
fprintf(fid,'%d\r\n',(constraintLabel(sum(modelparams.constraintDOF==[1 0 0 0 0 0],2)==6))');
fprintf(fid,'%s\r\n','dataend','fixboundary points XYT');
fprintf(fid,'%d\r\n',(constraintLabel(sum(modelparams.constraintDOF==[1 1 0 0 0 0],2)==6))');
fprintf(fid,'%s\r\n','dataend','fixboundary points XYZT');
fprintf(fid,'%d\r\n',(constraintLabel(sum(modelparams.constraintDOF==[1 1 1 0 0 0],2)==6))');
fprintf(fid,'%s\r\n','dataend','fixboundary points XYZTZR');
fprintf(fid,'%d\r\n',(constraintLabel(sum(modelparams.constraintDOF==[1 1 1 0 0 1],2)==6))');
fprintf(fid,'%s\r\n','dataend','fixboundary points ALL');
fprintf(fid,'%d\r\n',(constraintLabel(sum(modelparams.constraintDOF==[1 1 1 1 1 1],2)==6))');
fprintf(fid,'%s\r\n',batchEnd{:});
fclose(fid);


% Run model with no load to determine nodal positions
% Run ADINA to load batch file
!c:\adina93\x64\aui.exe -b -m 800mb ConicalSubstrate.in
pause(0.1);
!c:\adina93\x64\adina.exe -b -s -mm 1000mw -t 3 ConicalSubstrate.dat
pause(0.1);

% List coordinates
!c:\adina93\x64\aui.exe -b -m 800mb ConicalCoordData.plo
pause(0.1);

% Import data
headerLines = 11;
fid = fopen('ConicalCoordData.txt');
C = textscan(fid, '%*s %f %f %f %f','Delimiter',' ', 'MultipleDelimsAsOne',1, 'HeaderLines',headerLines); 
fclose(fid);
FEcoords = [C{1,1} C{1,2} C{1,3} C{1,4}];

% List deflections
!c:\adina93\x64\aui.exe -b -m 800mb ConicalDispData.plo
pause(0.1)

% Import results
nHeaderLines = 11;
% Columns in results file
% a-displacement: theta-displacement
% b-displacement: r-displacement
% c-displacement: z-displacement
% v1-rotation: drdtheta
% v2-rotation: drdz
fid = fopen('ConicalDispData.txt');
C = textscan(fid,'%*s %*d %f %f %f %f %f\r\n','delimiter',' ','multipledelimsasone',1,'headerlines',nHeaderLines);
fclose(fid);
FEdisp = [C{1,1} C{1,2} C{1,3} C{1,4} C{1,5}];


% Record coordinates and displacements for primary mirror
xFE = FEcoords(:,2);
yFE = FEcoords(:,3);
zFE = FEcoords(:,4);
[tFE,rFE] = cart2pol(xFE,yFE);
dxFE = FEdisp(:,1);
dyFE = FEdisp(:,2);
dzFE = FEdisp(:,3);
drdtFE = FEdisp(:,4);
drdzFE = FEdisp(:,5);

if modelparams.realignFlag==true
    % Remove rigid body displacement (rb) terms
    Arb = [ones(size(yFE)),yFE,zFE];
    dxflat = dxFE - Arb*pinv(Arb)*dxFE;
    % Calculate axial and aximuthal slopes
    drdz = drdzFE - mean(drdzFE);
    drdt = drdtFE - mean(drdtFE);
else
    dxflat = dxFE;
    drdz = drdzFE;
    drdt = drdtFE;
end

% Calculate radial displacement
dr = dxflat.*cos(tFE) + dyFE.*sin(tFE);

FEdata = struct('x',xFE,'y',yFE,'z',zFE,'t',tFE,'r',rFE,...
    'dx',dxFE,'dy',dyFE,'dz',dzFE,'dr',dr,...
    'drdz',drdzFE,'drdt',drdtFE);


% Remove the edges of the slope data before calculating the interpolants
dT = geoparams.Wt/(2*nDivTheta);
dZ = geoparams.Lp/(2*nDivZ);
drdzValid = tFE>=min(tFE)+dT/2 & tFE<=max(tFE)-dT/2;
drdtValid = zFE>=min(zFE)+dZ/2 & zFE<=max(zFE)-dZ/2;

% Interpolate the FE data onto a grid
[T0,Z0] = meshgrid(linspace(-geoparams.Wt/2,geoparams.Wt/2,2*nDivTheta+1),...
                   linspace(-geoparams.Lp/2,geoparams.Lp/2,2*nDivZ+1));
[~,idx] = sortrows([round(tFE/dT)*dT,round(zFE/dZ)*dZ]);
drgrid = reshape(dr(idx),[121,81]);
drdzgrid = reshape(drdz(idx),[121,81]);
drdtgrid = reshape(drdt(idx),[121,81]);
% drgrid = griddata(tFE,zFE,dr,T0,Z0);
% drdzgrid = griddata(tFE,zFE,drdz,T0,Z0);
% drdtgrid = griddata(tFE,zFE,drdt,T0,Z0);


% Interpolate deformation and slopes onto grid
% drFEinterp = scatteredInterpolant(tFE,zFE,dr,'linear','linear');
% drdzFEinterp = scatteredInterpolant(tFE(drdzValid),zFE(drdzValid),drdz(drdzValid),'linear','linear');
% drdtFEinterp = scatteredInterpolant(tFE(drdtValid),zFE(drdtValid),drdt(drdtValid),'linear','linear');
drFEinterp = griddedInterpolant(T0',Z0',drgrid','linear','linear');
drdzFEinterp = griddedInterpolant(T0(:,2:end-1)',Z0(:,2:end-1)',drdzgrid(:,2:end-1)','linear','linear');
drdtFEinterp = griddedInterpolant(T0(:,2:end-1)',Z0(:,2:end-1)',drdtgrid(:,2:end-1)','linear','linear');

FEinterp = struct('dr',drFEinterp,'drdz',drdzFEinterp,'drdt',drdtFEinterp);