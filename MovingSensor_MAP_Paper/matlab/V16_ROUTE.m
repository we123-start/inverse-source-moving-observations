function R=V16_ROUTE(truth,kind,noiseCase,noise,cfg)
%V16_ROUTE Paper Algorithm 2: measured-flux direction and recursive times.
% There are NO preassigned future observation windows. After each movement,
% the next 80-sample window begins on the first grid node not before arrival.
% The direction uses primary-endpoint and two auxiliary noisy flux values.
% Auxiliary readings never enter MAP reconstruction data.
% In the 12 paper experiments, MAP uses the accumulated observations after
% the route is generated. This is equivalent to a stagewise interleaving,
% because neither direction nor stopping decisions depend on the MAP fit.
% The fixed noise fixtures support four windows and three decisions.
maxWindows=cfg.maxWindows;
initialStart=cfg.initialStart;
kind=lower(char(kind));noiseCase=lower(char(noiseCase));
nt=cfg.mainPerWindow; dt=cfg.dt;
assert(maxWindows==floor(maxWindows) && maxWindows>=1 && maxWindows<=4, ...
    'The supplied fixed-noise fixtures permit at most four windows.');
assert(numel(noise.main)>=maxWindows*nt && ...
    numel(noise.nav)>=2*(maxWindows-1), ...
    'Insufficient primary or auxiliary noise samples.');
assert(all(isfinite(noise.main(:))) && all(isfinite(noise.nav(:))));
initialIndex=round(initialStart/dt);
assert(initialIndex>=1 && abs(initialIndex*dt-initialStart)<1e-12, ...
    'Initial sensing time must be a positive grid node.');
assert(initialIndex+nt-1<=cfg.steps, ...
    'Initial window extends beyond the forward solver time horizon.');
W=nan(maxWindows,2);
W(1,:)=[initialIndex,initialIndex+nt-1]*dt;
theta=nan(1,maxWindows);theta(1)=cfg.initial.(kind);
primaryEnd=zeros(1,maxWindows);primaryClean=zeros(1,maxWindows);
navMeasured=nan(maxWindows-1,2);navTimes=nan(1,maxWindows-1);
navCoords=nan(maxWindows-1,2);grad=nan(1,maxWindows-1);
dir=zeros(1,maxWindows-1); travel=zeros(1,maxWindows-1);
stepLength=zeros(1,maxWindows-1); reverse=false(1,maxWindows-1);
rawPeak=false(1,maxWindows-1);stoppedPeak=false(1,maxWindows-1);
arrival=nan(1,maxWindows-1);earliestNext=nan(1,maxWindows-1);
prev=0;nMoves=0;nNav=0;stopStage=0;stopReason='';
finishAfterReversal=false;
for k=1:maxWindows
    assert(all(isfinite(W(k,:))),'Missing recursively generated window %d.',k);
    a=theta(k); idx=round(W(k,1)/dt):round(W(k,2)/dt);
    assert(numel(idx)==nt && idx(end)<=cfg.steps, ...
        'Observation indices exceed the simulation horizon.');
    times=idx(:)*dt;noiseIdx=(k-1)*nt+(1:nt);
    clean=reviewer2_model('measuretruth',truth,a*ones(nt,1),times);
    measured=clean+noise.main(noiseIdx);
    primaryClean(k)=clean(end);primaryEnd(k)=measured(end);
    fprintf('%s stage %d: window [%d,%d] dt, position %.12g*pi/20\n', ...
        kind,k,idx(1),idx(end),a/(pi/20));
    if finishAfterReversal
        stopStage=k;stopReason='direction_reversal';break;
    end
    if k==maxWindows
        stopStage=k;stopReason='maximum_windows';break;
    end
    aux=mod(a+[-cfg.navOffset;cfg.navOffset],2*pi);
    zclean=reviewer2_model('measuretruth',truth,aux,times(end)*ones(2,1));
    zmeasured=zclean(:)+noise.nav(2*k-1:2*k);
    nNav=k;navCoords(k,:)=aux(:)';navTimes(k)=times(end);
    navMeasured(k,:)=zmeasured(:)';
    A=V16_MOVEMENT_DECISION(primaryEnd(k),zmeasured(1),zmeasured(2), ...
        prev,cfg.m.(kind),cfg.c1,cfg.rho0/2^(k-1), ...
        cfg.speed.(kind),cfg.navOffset,k>=cfg.minStopWindows);
    grad(k)=A.gradient;rawPeak(k)=A.threePointCandidate;
    stoppedPeak(k)=A.stopAtCurrent;
    if A.stopAtCurrent
        stopStage=k;stopReason='three_point_local_peak';break;
    end
    % Earliest possible next sensing window, after the actual boat travel.
    % No artificial reversal delay or preassigned future window is used.
    [sNext,eNext]=V16_NEXT_WINDOW(W(k,2),A.travelTime,dt,nt);
    if eNext>cfg.steps*dt+1e-12
        stopStage=k;stopReason='time_horizon';break;
    end
    theta(k+1)=theta(k)+A.direction*A.stepLength;
    W(k+1,:)=[sNext,eNext];
    dir(k)=A.direction;travel(k)=A.travelTime;
    stepLength(k)=A.stepLength;reverse(k)=A.reversed;
    arrival(k)=W(k,2)+travel(k);earliestNext(k)=sNext;
    nMoves=k;prev=A.direction;
    if A.reversed
        finishAfterReversal=true;
    end
end
assert(stopStage>=1 && stopStage<=maxWindows);
assert(nMoves==stopStage-1);
W=W(1:stopStage,:);theta=theta(1:stopStage);
% Construct the continuous, bounded-speed stop-and-go main-boat path.
pTime=[0;W(1,1);W(1,2)];
pTheta=theta(1)*ones(3,1);
for k=1:nMoves
    ts=[arrival(k);W(k+1,1);W(k+1,2)];
    as=theta(k+1)*ones(3,1);
    for j=1:3
        if ts(j)>pTime(end)+1e-14
            pTime(end+1,1)=ts(j); %#ok<AGROW>
            pTheta(end+1,1)=as(j); %#ok<AGROW>
        end
    end
end
if pTime(end)<cfg.steps*dt-1e-12
    pTime(end+1,1)=cfg.steps*dt;
    pTheta(end+1,1)=theta(end);
end
assert(all(diff(pTime)>0),'Nonincreasing sensor path times.');
assert(max(abs(diff(pTheta)./diff(pTime)))<=cfg.speed.(kind)+1e-7, ...
    'Sensor path exceeds the main boat speed limit.');
actual=round(theta/cfg.angleUnit);
expected=cfg.figIndex.(kind);
figureMatch=(stopStage==numel(expected) && isequal(actual,expected));
R=struct('kind',kind,'noiseCase',noiseCase,'mode','paper', ...
    'angles',theta,'anglesWrapped',mod(theta,2*pi), ...
    'displayIndex',actual,'expectedFigureIndex',expected, ...
    'matchesFigure',figureMatch,'mainWindows',W, ...
    'primaryMeasured',primaryEnd(1:stopStage), ...
    'primaryClean',primaryClean(1:stopStage), ...
    'navMeasured',navMeasured(1:nNav,:), ...
    'navCoords',navCoords(1:nNav,:), ...
    'navTimes',navTimes(1:nNav), ...
    'estimatedAngularDerivative',grad(1:nNav), ...
    'threePointCandidate',rawPeak(1:nNav), ...
    'threePointStop',stoppedPeak(1:nNav), ...
    'direction',dir(1:nMoves),'stepLength',stepLength(1:nMoves), ...
    'travelTime',travel(1:nMoves), ...
    'arrivalTime',arrival(1:nMoves), ...
    'earliestNextStart',earliestNext(1:nMoves), ...
    'reversalDetected',reverse(1:nMoves), ...
    'nCompletedWindows',stopStage,'nMoves',nMoves, ...
    'mainSamples',stopStage*nt,'auxiliarySamples',2*nNav, ...
    'stopStage',stopStage,'stopReason',stopReason, ...
    'pathTime',pTime,'pathTheta',pTheta, ...
    'feasible',true,'mainSpeed',cfg.speed.(kind), ...
    'rho0',cfg.rho0,'c1',cfg.c1);
fprintf('%s: %d windows, %d primary, %d auxiliary, reason=%s\n', ...
    kind,stopStage,R.mainSamples,R.auxiliarySamples,stopReason);
end
