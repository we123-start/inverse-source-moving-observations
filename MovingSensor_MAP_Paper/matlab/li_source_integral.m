function [F,dF] = li_source_integral(model,kind,z,b,nAng,nRad)
% Shape-fitted moving numerical quadrature for sharp characteristic chi_D.
% The original Q2 FEM stiffness, mass and backward Euler solver are reused.
% No smoothed source is used in the target posterior.
% dF is the derivative w.r.t. WHITENED parameter z.
if nargin<5 || isempty(nAng),nAng=96;end
if nargin<6 || isempty(nRad),nRad=8;end
z=z(:);d=numel(z);F=[];dF=[];
if strcmp(kind,'peanut')
    if nargout>1,[F,dF]=li_peanut_cutcell(model,z,b,16);
    else,F=li_peanut_cutcell(model,z,b,16);end
    return;
end
xi=reviewer2_shapes('physical',kind,z,'H1');
if ~reviewer2_shapes('valid',kind,xi,1e-6),return;end
[p,rho,qw]=local_quadrature(nAng,nRad);
N=numel(p);dx=zeros(d,N);dy=dx;dj=dx;
r0=xi(1);phi0=xi(2);s=xi(3);
dr0=1/(pi*(1+z(1)^2));dphi0=2/(1+z(2)^2);
ds=1/(pi*(1+z(3)^2));
switch lower(kind)
  case 'circle'
    X=cos(p);Y=sin(p);cross=ones(1,N);
  case 'fourleaf'
    a=.4+.2*cos(4*p);X=a.*cos(p);Y=a.*sin(p);cross=a.^2;
  case 'kite'
    X=cos(p)+.65*cos(2*p)-.65;Y=1.5*sin(p);
    Xp=-sin(p)-1.3*sin(2*p);Yp=1.5*cos(p);
    cross=X.*Yp-Y.*Xp;
  otherwise
    error('Unknown source shape');
end
x=r0*cos(phi0)+rho.*s.*X;
y=r0*sin(phi0)+rho.*s.*Y;
jac=rho*s*s.*cross;
dx(1,:)=dr0*cos(phi0);dy(1,:)=dr0*sin(phi0);
dx(2,:)=-r0*sin(phi0)*dphi0;dy(2,:)=r0*cos(phi0)*dphi0;
dx(3,:)=rho.*X*ds;dy(3,:)=rho.*Y*ds;
dj(3,:)=2*jac*ds/s;
r=hypot(x,y);
if any(~isfinite(r)) || any(r>=1-1e-13),return;end
ang=mod(atan2(y,x),2*pi);
hx=1/(2*model.m);hy=pi/model.n;
ri=min(model.m-1,floor(r/(2*hx)));
ai=min(model.n-1,floor(ang/(2*hy)));
u=(r-(2*ri+1)*hx)/hx;
v=(ang-(2*ai+1)*hy)/hy;
Lr=[u.*(u-1)/2;1-u.^2;u.*(u+1)/2];
La=[v.*(v-1)/2;1-v.^2;v.*(v+1)/2];
Dr=[u-.5;-2*u;u+.5]/hx;
Da=[v-.5;-2*v;v+.5]/hy;
if nargout>=2
    rr=max(r,1e-12);
    rder=(dx.*repmat(x,d,1)+dy.*repmat(y,d,1))./repmat(rr,d,1);
    ader=(dy.*repmat(x,d,1)-dx.*repmat(y,d,1))./repmat(max(r.^2,1e-24),d,1);
end
F=zeros(model.nFree,1);
if nargout>=2,dF=zeros(model.nFree,d);end
w=qw.*jac;
nr=2*model.m+1;
for da=0:2
    for dr=0:2
        radial=2*ri+dr;
        angular=mod(2*ai+da,model.na);
        indices=1+angular*(nr-2)+radial;
        indices(radial==0)=1;
        alive=(radial~=nr-1);
        ind=indices(alive);
        base=La(da+1,:).*Lr(dr+1,:);
        F=F+accumarray(ind(:),b*(base(alive).*w(alive)).',[model.nFree,1],@sum,0);
        if nargout>=2
            br=La(da+1,:).*Dr(dr+1,:);
            ba=Da(da+1,:).*Lr(dr+1,:);
            for j=1:d
                derivative=(br.*rder(j,:)+ba.*ader(j,:)).*w+base.*qw.*dj(j,:);
                dF(:,j)=dF(:,j)+accumarray(ind(:),b*derivative(alive).',[model.nFree,1],@sum,0);
            end
        end
    end
end
end

function [angle,rho,w]=local_quadrature(nAng,nRad)
persistent cachedA cachedR cachedAngle cachedRho cachedW
if isequal(cachedA,nAng)&&isequal(cachedR,nRad)
    angle=cachedAngle;rho=cachedRho;w=cachedW;return;
end
[p,pa]=local_gauss(6);[r,rw]=local_gauss(nRad);
nb=ceil(nAng/6);ang=[];aw=[];
for k=0:nb-1
    ang=[ang,((k+.5)+p.'/2)*(2*pi/nb)]; %#ok<AGROW>
    aw=[aw,pa.'*(pi/nb)]; %#ok<AGROW>
end
ang=ang(:)';aw=aw(:)';
r=(r+1)/2;rw=rw/2;
angle=kron(ang,ones(1,nRad));rho=repmat(r.',1,numel(ang));
w=kron(aw,rw.');
cachedA=nAng;cachedR=nRad;cachedAngle=angle;cachedRho=rho;cachedW=w;
end

function [x,w]=local_gauss(n)
k=(1:n-1)';off=k./sqrt(4*k.^2-1);
[V,D]=eig(diag(off,1)+diag(off,-1));
[x,perm]=sort(diag(D));w=2*(V(1,perm)').^2;
x=x(:);w=w(:);
end
