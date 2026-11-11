function cfg=V16_CONFIG()
% Single production configuration for Q2/BE + MAP/SQP paper experiments.
% All observation windows arise recursively from the flux-driven route.
% The observation positions are j*pi/20 (not j*pi/40).
cfg=struct();
cfg.dt=1/400;cfg.steps=400;cfg.truthMesh=23;cfg.invMesh=20;
cfg.maxWindows=4;cfg.initialStart=cfg.dt;
cfg.figIndex=struct('circle',[26 16 6 11], ...
    'kite',[26 16 6 11], 'fourleaf',[29 19 9 14], ...
    'peanut',[4 19 34 27]);
cfg.angleUnit=pi/20;
cfg.initial=struct('circle',26*pi/20,'kite',26*pi/20, ...
    'fourleaf',29*pi/20,'peanut',4*pi/20);
cfg.m=struct('circle',10,'kite',10,'fourleaf',10,'peanut',15);
cfg.speed=struct('circle',20*pi,'kite',20*pi,'fourleaf',20*pi,'peanut',30*pi);
cfg.c1=1/20;cfg.rho0=1e-10;cfg.navOffset=pi/20;
cfg.sigma=struct('circle',0.05,'kite',0.05,'fourleaf',0.05,'peanut',0.01);
% The algorithm uses two auxiliary FLUX readings per decision, but no
% auxiliary-vessel kinematics or speed is included in the numerical model.
% Only the main-sensor path is constrained by cfg.speed.(kind).
% Require two distinct stationary observation sites before early stopping.
% This prevents a one-site stop conflicting with Theorem 2.1.
cfg.minStopWindows=2;
cfg.mainPerWindow=80;cfg.nAngInv=96;cfg.nRadInv=8;
cfg.baseSeed=20261101;
% Single, fixed initialization in the Gaussian/reference coordinates w.
% This is NOT the physical parameter xi. One shared 3D starting vector for
% circle, kite, fourleaf; one 5D vector for peanut. Chosen a priori, not
% inspecting MAP reconstruction errors.
cfg.mapW0Three=[0;0;-2];
% Strictly admissible noncircular five-parameter starting shape.
% In physical coordinates xi=[.6;.1;.1;.05;.05] (H1 scaling).
cfg.mapW0Peanut=[0.6;0.1;0.1;0.1;0.1];
cfg.mapInitial=struct('circle',cfg.mapW0Three,'kite',cfg.mapW0Three,...
    'fourleaf',cfg.mapW0Three,'peanut',cfg.mapW0Peanut);
end
