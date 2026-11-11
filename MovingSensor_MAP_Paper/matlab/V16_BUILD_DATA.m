function [D,info]=V16_BUILD_DATA(truth,kind,R,N,cfg)
% Recompute all seven comparison designs on the SAME 320 time instants.
% Truth-informed A/B/C are explicit synthetic stress-test controls.
% No use of reconstruction error or MAP results in sensor selection.
assert(numel(N.main)==320);
assert(R.nCompletedWindows==cfg.maxWindows, ...
    'Paper comparisons require four completed observation windows.');
W=R.mainWindows;times=zeros(cfg.mainPerWindow*cfg.maxWindows,1);
for k=1:4
    ii=(k-1)*cfg.mainPerWindow+(1:cfg.mainPerWindow);
    inds=round(W(k,1)/cfg.dt):round(W(k,2)/cfg.dt);
    assert(numel(inds)==cfg.mainPerWindow, ...
        'Each observation window needs 80 BE nodes');
    times(ii)=inds(:)*cfg.dt;
end
% Compare average |exact flux| for 160 candidate fixed sites, sampled at
% 80 of the common 320 times. This stage is truth-INFORMED by design.
candidate=(0:159)'*(2*pi/160);
trainingTimes=times(1:4:end); % exactly 80 of 320 time indices
score=zeros(160,1);
for j=1:160
    aa=candidate(j)*ones(numel(trainingTimes),1);
    f=reviewer2_model('measuretruth',truth,aa,trainingTimes);
    score(j)=mean(abs(f));
end
[~,order]=sort(score,'ascend');
low=[];lowScore=[];
for j=1:numel(order)
    a=candidate(order(j));
    if isempty(low) || all(abs(atan2(sin(a-low),cos(a-low)))>=pi/3-1e-12)
        low=[low a]; %#ok<AGROW>
        lowScore=[lowScore score(order(j))]; %#ok<AGROW>
    end
    if numel(low)==3,break;end
end
assert(numel(low)==3,'Not enough angularly separated low-flux controls');
% Nonadaptive two-site trajectory. Fixed BEFORE optimization.
offset=pi*[1/12+sqrt(2)*1e-5, ...
          -(17/180+sqrt(2)*1e-5), ...
          19/180+sqrt(2)*1e-5];
assert(all(abs(offset)/cfg.speed.(kind)<W(3,1)-W(2,2)));
policies={'adaptive','fixed_initial','fixed_low_A','fixed_low_B',...
          'local_A','local_B','local_C'};
D=struct();
for j=1:numel(policies)
    policy=policies{j};
    switch policy
        case 'adaptive',pos=R.angles;
        case 'fixed_initial',pos=cfg.initial.(kind)*ones(1,4);
        case 'fixed_low_A',pos=low(1)*ones(1,4);
        case 'fixed_low_B',pos=low(2)*ones(1,4);
        otherwise
            k=j-4;pos=[low(k),low(k),low(k)+offset(k),low(k)+offset(k)];
    end
    a=repelem(mod(pos(:),2*pi),cfg.mainPerWindow);
    clean=reviewer2_model('measuretruth',truth,a,times);
    nAuxPolicy=0;
    if strcmp(policy,'adaptive'), nAuxPolicy=R.auxiliarySamples; end
    D.(policy)=struct('kind',kind,'policy',policy,'noiseCase',R.noiseCase,'mode','algorithm',...
        'angles',a,'times',times,'positions',pos,'clean',clean,...
        'noise',N.main(:),'observed',clean+N.main(:),...
        'nMainObservations',320,'nAuxNavigationObservations',nAuxPolicy,...
        'nav',R,'rep',N.rep,'noiseSeed',N.seed);
end
for j=2:7
    assert(isequal(D.(policies{j}).times,D.adaptive.times));
    assert(isequal(D.(policies{j}).noise,D.adaptive.noise));
end
info=struct('kind',kind,'policies',{policies},'noiseCase',R.noiseCase,'mode','algorithm',...
     'lowFluxAngles',low,'lowFluxScores',lowScore,'gridAngles',candidate,...
     'gridFluxScores',score,'fixedSelectionTimes',trainingTimes,...
     'localOffsets',offset,'mainObservationTimes',times,...
     'mainObservationsPerPolicy',320,'navObservationsAdaptive',6,...
     'noiseOrigin',N.origin,'noiseSeed',N.seed,...
     'selectionTruthInformed',true);
end
