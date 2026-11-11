function [val,g,pred] = li_phi_grad(z,post)
% Returns negative log likelihood Phi, its gradient and model predictions.
[val,g,pred]=deal(inf,zeros(numel(z),1),[]);
if nargout>=2
    [f,df]=li_source_integral(post.model,post.kind,z,post.b,post.nAng,post.nRad);
else
    f=li_source_integral(post.model,post.kind,z,post.b,post.nAng,post.nRad);
end
if isempty(f),return;end
pred=post.H*f;
res=pred-post.y;
val=(res.'*res)/(2*post.sigma^2);
if nargout>=2,g=df.'*(post.H.'*res)/post.sigma^2;end
end
