function A=V16_MOVEMENT_DECISION(primaryEnd,minusFlux,plusFlux,prev,m,c1,rho,speed,navOffset,allowStop)
%V16_MOVEMENT_DECISION Pure decision logic for manuscript Algorithm 2.
% First test the three-point stopping criterion. Move ONLY if it fails.
% The auxiliary readings are not added to the inverse data likelihood.
% Inputs are measured (i.e. noisy when appropriate) signed fluxes.
% allowStop must be false until TWO distinct stationary sites have been
% observed, consistent with the theorem's two-location requirement.
if nargin<10,allowStop=true;end
assert(isscalar(allowStop) && (islogical(allowStop) || isnumeric(allowStop)));
assert(all(isfinite([primaryEnd,minusFlux,plusFlux,c1,rho,speed,navOffset])));
assert(isscalar(m) && m>=1 && m==floor(m));
assert(c1>0 && rho>0 && speed>0 && navOffset>0);
A=struct('primaryFlux',primaryEnd,'minusFlux',minusFlux,'plusFlux',plusFlux, ...
    'primaryMagnitude',abs(primaryEnd),'auxiliaryMax',max(abs(minusFlux),abs(plusFlux)), ...
    'gradient',(abs(plusFlux)-abs(minusFlux))/(2*navOffset), ...
    'threePointCandidate',false,'eligibleStop',logical(allowStop), ...
    'stopAtCurrent',false,'direction',0,'reversed',false, ...
    'stepLength',0,'travelTime',0,'rhoApplied',0);
A.threePointCandidate=A.primaryMagnitude>A.auxiliaryMax;
A.stopAtCurrent=logical(allowStop) && A.threePointCandidate;
if A.stopAtCurrent
    return;
end
% Positive angular derivative -> counterclockwise. Zero -> clockwise.
A.direction=2*(A.gradient>0)-1;
A.reversed=(prev~=0 && A.direction~=prev);
if prev==0 || ~A.reversed
    nominal=m*c1*pi;
else
    nominal=floor(m/2)*c1*pi;
end
A.rhoApplied=rho;
A.stepLength=nominal+rho;
A.travelTime=A.stepLength/speed;
end
