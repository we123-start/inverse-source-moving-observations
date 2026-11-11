function [F,dF] = li_peanut_cutcell(model,z,b,thetaOrder)
% VECTORIZED radial-exact Q2 source load for the sharp Fourier peanut.
% Same FEM stiffness/mass and BE flux transfer as the submitted code.
% All candidate q(theta) are integrated radially exactly over every Q2
% radial cell, while angular integrals use Gauss-Legendre thetaOrder nodes
% per FEM angular cell. No likelihood smoothing or fixed-point jump grid.
% Gradient is analytic almost everywhere in parameter space.
%
% For M=2: xi=[z1,z2,z3,z4/2,z5/2] (Gaussian prior B as in paper),
% q(theta)=xi1/2+xi2*cos(theta)+xi3*sin(theta)
%                 +xi4*cos(2theta)+xi5*sin(2theta).
if nargin<4,thetaOrder=16;end
persistent cache
key=[model.m model.n model.nFree thetaOrder];
if isempty(cache) || ~isequal(cache.key,key)
    cache=build_cache(model,thetaOrder);
end
z=z(:);d=numel(z);assert(d==5);
xi=reviewer2_shapes('physical','peanut',z,'H1');
if ~reviewer2_shapes('valid','peanut',xi,1e-6)
    F=[];if nargout>1,dF=[];end;return;
end
m=model.m;n=model.n;nt=thetaOrder;hx=1/(2*m);
% Columns of coeffRad: angular element a, radial element r, angular
% Gauss point t, matching Q columns 3*((a*m+r)*nt+t)+dr.
coeffRad=cache.coeffRad;                 % [5 x (n*m*nt)]
q=z.'*coeffRad;                         % [1 x K]
sRaw=(q-cache.midRad)/hx;
s=max(-1,min(1,sRaw));
s2=s.*s;s3=s2.*s;s4=s2.*s2;
Ir=[s3/6-s2/4+5/12; s-s3/3+2/3; s3/6+s2/4-1/12];
Kr=[s4/8-s3/6-7/24; s2/2-s4/4-1/4; s4/8+s3/6+1/24];
radialIntegrals=hx*(bsxfun(@times,Ir,cache.midRad)+hx*Kr);
F=b*(cache.Q*radialIntegrals(:));
if nargout>1
    L=[s.*(s-1)/2;1-s2;s.*(s+1)/2];
    % d/dq int_0^q N_i(r) r dr = q*N_i(q) inside this cell.
    active=(sRaw>-1 & sRaw<1);
    jac=(cache.midRad+hx*s).*active;
    local=bsxfun(@times,L,jac);
    dh=zeros(3*numel(q),d);
    for j=1:d
        dh(:,j)=reshape(bsxfun(@times,local,coeffRad(j,:)),[],1);
    end
    dF=b*(cache.Q*dh);
end
end

function S=build_cache(model,nt)
m=model.m;n=model.n;hx=1/(2*m);hy=pi/n;
[g,w]=local_gauss(nt);
nr=2*m+1;
% 5 x nt x n, then duplicate per radial cell.
C=zeros(5,nt,n);
for a=0:n-1
    theta=(2*a+1)*hy+hy*g';
    C(:,:,a+1)=[.5*ones(1,nt);cos(theta);sin(theta);...
                 .5*cos(2*theta);.5*sin(2*theta)];
end
S.coeffRad=reshape(repmat(reshape(C,[5 nt 1 n]),[1 1 m 1]),5,[]);
mid=(2*(0:m-1)+1)*hx;
S.midRad=repmat(repelem(mid,nt),1,n);
% Fixed radial-Q2/angular-Q2 assembly map Q: Q*h integrates one source.
N=3*n*m*nt;
rows=zeros(3*N,1);cols=rows;val=rows;used=0;
for a=0:n-1
    for r=0:m-1
        for t=1:nt
            v=g(t);La=[v*(v-1)/2;1-v^2;v*(v+1)/2];
            for dr=0:2
                col=3*((a*m+r)*nt+(t-1))+dr+1;
                ir=2*r+dr;
                if ir==nr-1,continue;end % homogeneous r=1 boundary
                for da=0:2
                    if ir==0,id=1;
                    else,id=1+mod(2*a+da,2*n)*(nr-2)+ir;
                    end
                    used=used+1;
                    rows(used)=id;cols(used)=col;
                    val(used)=hy*w(t)*La(da+1);
                end
            end
        end
    end
end
S.Q=sparse(rows(1:used),cols(1:used),val(1:used),model.nFree,N);
S.key=[m n model.nFree nt];
end
function [x,w]=local_gauss(n)
k=(1:n-1)';off=k./sqrt(4*k.^2-1);
[V,D]=eig(diag(off,1)+diag(off,-1));
[x,ix]=sort(diag(D));w=2*(V(1,ix)').^2;
x=x(:);w=w(:);
end
