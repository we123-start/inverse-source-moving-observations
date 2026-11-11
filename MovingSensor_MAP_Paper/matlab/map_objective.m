function [f,g]=map_objective(w,post)
w=w(:);[phi,gphi]=li_phi_grad(w,post);
if isfinite(phi)
    f=phi+0.5*(w'*w);g=gphi(:)+w;
else
    [c,~]=map_constraints(w,post.kind);
    violation=sum(max(0,c(:)).^2);
    f=1e7+1e7*violation+0.5*(w'*w);
    g=w;
    if violation>0
        % SQP was allowed to evaluate outside the physical domain. Match
        % the gradient of this finite extension to its objective value.
        for jj=1:numel(w)
            step=1e-5*max(1,abs(w(jj)));
            wp=w;wm=w;wp(jj)=wp(jj)+step;wm(jj)=wm(jj)-step;
            cp=map_constraints(wp,post.kind);
            cm=map_constraints(wm,post.kind);
            vp=sum(max(0,cp(:)).^2);
            vm=sum(max(0,cm(:)).^2);
            g(jj)=g(jj)+1e7*(vp-vm)/(2*step);
        end
    end
end
end
