function d = MAP_STAGE_DATA(dFull,k,cfg)
%MAP_STAGE_DATA Only the first k stationary observation windows.
assert(isscalar(k) && k >= 1 && k <= 4 && k == floor(k));
n = k*cfg.mainPerWindow;
assert(numel(dFull.observed) == 4*cfg.mainPerWindow);
d = dFull;
d.angles = dFull.angles(1:n);
d.times = dFull.times(1:n);
d.clean = dFull.clean(1:n);
d.noise = dFull.noise(1:n);
d.observed = dFull.observed(1:n);
d.positions = dFull.positions(1:k);
d.nMainObservations = n;
end
