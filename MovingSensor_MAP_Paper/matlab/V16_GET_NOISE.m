function N=V16_GET_NOISE(root,kind,noiseCase,rep,cfg)
% Exact Python/MATLAB paired Gaussian values in first noisy repetition.
if nargin<5,cfg=V16_CONFIG();end
assert(rep>=1 && rep==round(rep));
noiseCase=lower(char(noiseCase));
assert(any(strcmp(noiseCase,{'zero','noisy'})));
N=struct('nav',zeros(6,1),'main',zeros(320,1),...
    'noiseCase',noiseCase,'rep',rep,'seed',0,'origin','deterministic zero noise');
if strcmp(noiseCase,'zero'),return;end
sigma=cfg.sigma.(kind);
if rep==1
    f=fullfile(root,'fixtures',['V16_REPRO_NOISE_',kind,'.mat']);
    assert(exist(f,'file')==2,'Missing exact cross-language noise: %s',f);
    % Do not import metadata named 'case': it is a MATLAB keyword.
    % Loading only these arrays also avoids a harmless keyword warning.
    B=load(f,'navZ','mainZ','navSeed');
    assert(numel(B.navZ)==6 && numel(B.mainZ)==320,'Invalid noise length');
    N.nav=sigma*B.navZ(:);N.main=sigma*B.mainZ(:);
    N.seed=double(B.navSeed);
    N.origin='fixed NumPy PCG64 draws stored as MAT; identical in MATLAB';
else
    shapes={'circle','kite','fourleaf','peanut'};id=find(strcmp(kind,shapes),1);
    seed=cfg.baseSeed+100000*rep+1000*id;
    prev=rng;rng(seed,'twister');
    N.nav=sigma*randn(6,1);N.main=sigma*randn(320,1);
    rng(prev);
    N.seed=seed;N.origin='independent MATLAB twister replication';
end
assert(all(isfinite(N.nav)) && all(isfinite(N.main)));
end
