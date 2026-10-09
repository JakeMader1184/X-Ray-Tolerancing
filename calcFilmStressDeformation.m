% function calcFilmStressDeformation(geoparams, constraintVect, Nfront, Nback)
% Calculates the deformation of a mirror resulting from film stress applied
% to the surface.
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
% constraints - a matrix (size Nx5) containing the positions of each
%               constraint, and which degrees of freedom are fixed.
%               [tC, zC, DOFxC, DOFyC, DOFzC]:
%               tC, zC - azimuthal [rad] and z [m] coordinates of
%                        constraint point.
%               DOFxC, DOFyC, DOFzC - flag for whether degree of freedom is
%                        fixed (1) or free (0).
%
% OUTPUTS:
addpath('..\Ray trace functions')
clear all
clc


hs = 0.5e-3;
Lp = 100e-3; % Length of primary mirror
Wmax = Lp; % Maximum mirror arc-width.

hf = 0.02e-6;

% r0Vect = [250:250:1500]*1e-3;
r0Vect = 1500e-3;
% Define stress magnitudes [N/m]. 9 different stress distributions will
% be tested, defined by Legendre-Legendre polynomials, which each have
% a peak value defined by this parameter.
stressFront = -0.2;

t0 = 0*pi/180; % The mirror must not cross +-pi, because this is where there is a discontinuity in the atan2 function
z0 = 10.0;

% Set number of elements for the FE model
nDivTheta = 40;             % Number of elements in the azimuthal direction
nDivZ = 60;                 % Number of elements in the axial direction

% Set substrate material parameters
E = 169e9;                  % Elastic modulus [Pa]
nu = 0.064;                 % Poisson's ratio

% Define constraint positions and degrees of freedom
constraintPos = [0,0;
    0,Lp/2];
constraintDOF = [1,1,1,0,0,1;
                 1,1,0,0,0,0];
realignFlag=true;

% dzpost = 0e-3;
% constraintPos = [Lp/4/r0,Lp/4+dzpost;
%                  -Lp/4/r0,Lp/4+dzpost;
%                  -Lp/4/r0,-Lp/4+dzpost;
%                  Lp/4/r0,-Lp/4+dzpost];
% constraintDOF = [1,1,1,0,0,0;
%                  1,1,1,0,0,0;
%                  1,1,1,0,0,0;
%                  1,1,1,0,0,0];
% realignFlag=false;



for idx = 1:length(r0Vect)
    r0 = r0Vect(idx);
    Wt = 2*pi/floor(2*pi*r0/Wmax/2)/2; % Angular width [rad]
    
    % Create input structures
    geoparams = struct('r0',r0,'t0',t0,'z0',z0,'Lp',Lp,'Wt',Wt,'zgap',0,...
        'Npointsz',2*nDivZ+1,'Npointst',2*nDivTheta+1,...
        'eP',[0 0 0],'eS',[0 0 0],'phiP',[0 0 0],'phiS',[0 0 0]);
    matparams = struct('E',E,'nu',nu);
    modelparams = struct('nDivTheta',nDivTheta,'nDivZ',nDivZ,...
        'hs',hs,'hf',hf,'hb',0,...
        'stressFront',stressFront,...
        'constraintPos',constraintPos,...
        'constraintDOF',constraintDOF,'realignFlag',realignFlag);
    tic
    mirrorPosition = 'P';
    % Run finite element model for primary mirror
    [FEdataP, FEinterpP] = calcConicalDeformation(geoparams, matparams, modelparams,mirrorPosition);
    
    toc
    
    tic
    mirrorPosition = 'S';
    % Run finite element model for secondary mirror
    [FEdataS, FEinterpS] = calcConicalDeformation(geoparams, matparams, modelparams,mirrorPosition);
    
    toc
    
    mirrorName = ['T' num2str(round(hs*1e6)) 'R' num2str(r0*1000,'%04d') 'L' num2str(round(Lp*1e3)) 'Nf' num2str(abs(stressFront)*10,'%02d')];
    
    save(['Deformation data\' mirrorName '.mat'],...
        'geoparams','matparams','modelparams','FEdataP','FEdataS', 'FEinterpP', 'FEinterpS');
end