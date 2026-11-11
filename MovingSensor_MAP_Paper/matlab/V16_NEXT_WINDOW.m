function [nextStart,nextEnd,startIndex,endIndex]=V16_NEXT_WINDOW(currentEnd,travelTime,dt,nSamples)
%V16_NEXT_WINDOW Recursive grid-aligned sensing window after a sensor move.
% No future sensing time is prescribed. For an 80-sample window, the new
% interval has endpoints [nextStart, nextStart+79*dt].
assert(isscalar(currentEnd) && isfinite(currentEnd) && currentEnd>=0);
assert(isscalar(travelTime) && isfinite(travelTime) && travelTime>0);
assert(isscalar(dt) && isfinite(dt) && dt>0);
assert(isscalar(nSamples) && nSamples==floor(nSamples) && nSamples>=1);
arrival=currentEnd+travelTime;
startIndex=ceil(arrival/dt);
% Safeguard against a roundoff-induced start just before actual arrival.
while startIndex*dt < arrival
    startIndex=startIndex+1;
end
endIndex=startIndex+nSamples-1;
nextStart=startIndex*dt;
nextEnd=endIndex*dt;
assert(nextStart>=arrival);
end
