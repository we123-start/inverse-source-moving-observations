function R=MAP_SOLVE_MEMORY(invModel,kind,D,w0,cfg)
% Compute one MAP solution entirely in memory; NO MAT files or figures saved.
% Solve the ORIGINAL Q2/BE inverse problem, not a regression surrogate.
% Fixed one-start SQP in reference coordinate w; final refinement retained.
fullClock=tic; % Whole inverse-case computation (not sensor travel time)
Hfull=reviewer2_model('design',invModel,D.angles,D.times);
[U,S,~]=svd(Hfull,'econ');ss=diag(S);
r=max(1,sum(ss>1e-12*ss(1)));
post=struct('model',invModel,'kind',kind,'H',U(:,1:r)'*Hfull, ...
    'y',U(:,1:r)'*D.observed,'b',reviewer2_shapes('amplitude',kind), ...
    'sigma',cfg.sigma.(kind),'nAng',cfg.nAngInv,'nRad',cfg.nRadInv, ...
    'rank',r,'angles',D.angles,'times',D.times);
optClock=tic;
[w,obj,trials]=solve_map_sqp(post,w0);
optWallSeconds=toc(optClock);
xi=reviewer2_shapes('physical',kind,w,'H1');
assert(reviewer2_shapes('valid',kind,xi,1e-6), ...
    'Final MAP violates shape admissibility');
% Publish the FULL data-space objective, not the projected objective.
% The projection is only a low-rank numerical optimization accelerator;
% its removed orthogonal noise term depends on the sensor policy.
F=li_source_integral(invModel,kind,w,reviewer2_shapes('amplitude',kind), ...
    cfg.nAngInv,cfg.nRadInv);
assert(~isempty(F),'Invalid MAP source integral.');
fullResidual=Hfull*F-D.observed;
fullObj=.5*sum((fullResidual/cfg.sigma.(kind)).^2)+.5*(w'*w);
truth=reviewer2_shapes('truth',kind);
% Angle is circular, so use shortest angular error in the parameter norm.
d=xi(:)-truth(:);
if ~strcmp(kind,'peanut'),d(2)=atan2(sin(d(2)),cos(d(2)));end
relPar=norm(d)/norm(truth);
% Unit-disk grid approximation of relative symmetric difference and IoU.
q=linspace(-1,1,321);[X,Y]=meshgrid(q,q);
inside=X.^2+Y.^2<1;
P=logical(reviewer2_shapes('mask',kind,truth,X(:),Y(:)));
Q=logical(reviewer2_shapes('mask',kind,xi,X(:),Y(:)));
P=P & inside(:);Q=Q & inside(:);
relSym=sum(xor(P,Q))/max(1,sum(P));iou=sum(P&Q)/max(1,sum(P|Q));
R=struct('kind',kind,'policy',D.policy,'mode',D.mode,'rep',D.rep, ...
    'wMAP',w,'xiMAP',xi,'initialW',w0(:),'objective',fullObj, ...
    'optimizationObjectiveProjected',obj,'trials',{trials}, ...
    'relativeParameterError',relPar, ...
    'relativeSymmetricDifference',relSym,'IoU',iou, ...
    'nMainObservations',numel(D.observed),'nStarts',1, ...
    'SVDprojectionRank',r,'positions',D.positions,'times',D.times, ...
    'noiseSeed',D.noiseSeed, ...
    'MAPoptimizerElapsedSeconds',optWallSeconds, ...
    'MAPtotalWallSeconds',toc(fullClock));
end
