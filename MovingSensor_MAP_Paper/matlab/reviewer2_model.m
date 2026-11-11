function varargout=reviewer2_model(action,varargin)
% Only run the path check once. RUN_MAIN adds all paths automatically.
persistent femPathReady
if isempty(femPathReady)
    originalFemDir=fullfile(fileparts(mfilename('fullpath')),'original_fem');
    assert(exist(originalFemDir,'dir')==7, ...
        'Missing original Q2 FEM helpers: %s',originalFemDir);
    addpath(originalFemDir,'-begin');
    assert(exist('gaussian_nodes_weights_rect_ref','file')==2, ...
        'Missing gaussian_nodes_weights_rect_ref.m');
    femPathReady=true;
end
% Exact full-order BE input-output transfer from the user's Q2 polar FEM.
% The original mesh, quadrature, basis and stiffness assembly are preserved.
% The precomputed discrete BE response kernel avoids repeat time stepping;
% no eigenmode truncation or surrogate forward model is used.
switch lower(action)
    case 'setup'
        varargout{1}=setup_model(varargin{:});
    case 'truth'
        varargout{1}=build_truth(varargin{:});
    case 'measuretruth'
        varargout{1}=evaluate_truth(varargin{:});
    case 'design'
        varargout{1}=create_design(varargin{:});
    otherwise, error('Unknown reviewer2_model action');
end
end

function model=setup_model(m,n,dt,nStep,buildKernel,scheme,sourceOrder)
if nargin<5,buildKernel=true;end
if nargin<6 || isempty(scheme),scheme='BE';end
if nargin<7 || isempty(sourceOrder),sourceOrder=6;end
scheme=upper(char(scheme));
assert(strcmp(scheme,'BE'),'Reproducibility package supports backward Euler only');
fprintf('  Constructing original polar Q2 FEM: m=%d n=%d ...\n',m,n);
tic;
global gx gy gw
[gx,gy,gw]=gaussian_nodes_weights_rect_ref(6);
hx=1/(2*m);hy=pi/n;
[xy,elems,~]=supermeshstructurepolar(m,n);
[Aa,Atheta,Mraw]=superstiffmatrixpolar(xy,elems,hx,hy);
Ar=Aa+Atheta;
na=2*n;nr=2*m+1;nraw=size(xy,2); nFree=1+na*(nr-2);
% P identifies angular seam and ALL r=0 nodes; no artificial Dirichlet at
% the pole, and radial r=1 nodes are excluded (Dirichlet only there).
map=zeros(nraw,1);
for ia=0:na
    aa=mod(ia,na);
    for ir=0:(nr-2)
        raw=ia*nr+ir+1;
        if ir==0,map(raw)=1;
        else,map(raw)=1+aa*(nr-2)+ir;
        end
    end
end
sel=find(map>0);
P=sparse(sel,map(sel),1,nraw,nFree);
A=P'*Ar*P;M=P'*Mraw*P;
A=(A+A')/2;M=(M+M')/2;
assert(all(diag(M)>0),'FE mass matrix has nonpositive diagonal');
% True outward normal derivative at r=1: Q2 radial trace derivative.
I=zeros(na*2,1);J=I;V=I;p=0;
for ia=0:na-1
    for item=1:2
        if item==1,ir=nr-3; val=0.5/hx;
        else,ir=nr-2;val=-2/hx;
        end
        raw=ia*nr+ir+1;
        p=p+1; I(p)=ia+1;J(p)=map(raw);V(p)=val;
    end
end
H=sparse(I,J,V,na,nFree);
% Fast source load operator W: OPTIONAL higher source quadrature,
% without changing the mass, stiffness or backward-Euler scheme.
[qgx,qgy,qgw]=source_rectangle_gauss(sourceOrder);
phi=superbasicu(qgx',qgy',hx,hy);
nQuadPerEle=numel(qgx);nElems=size(elems,2);
nQuad=nQuadPerEle*nElems;
Iq=zeros(9*nQuad,1);Jq=Iq;Vq=Iq;qX=zeros(nQuad,1);qY=qX;
pos=0;
for e=1:nElems
    node=elems(:,e);
    r0=mean(xy(1,node));a0=mean(xy(2,node));
    rr=r0+hx*qgx; aa=a0+hy*qgy;
    inds=(e-1)*nQuadPerEle+(1:nQuadPerEle);
    qX(inds)=rr.*cos(aa);qY(inds)=rr.*sin(aa);
    wt=(qgw(:)'.*rr')*(hx*hy);
    for loc=1:9
        k=map(node(loc));
        if k==0,continue;end
        off=pos+(1:nQuadPerEle);
        Iq(off)=k;
        Jq(off)=inds;
        Vq(off)=phi(loc,:).*wt;
        pos=pos+nQuadPerEle;
    end
end
W=sparse(Iq(1:pos),Jq(1:pos),Vq(1:pos),nFree,nQuad);
% Backward Euler: (M+dt*A)U^{n+1}=M U^n+dt F.
L=M+dt*A;R=M;
% Build factor once: MATLAB base functionality. decomposition is available
% in modern MATLAB; fallback uses LU if unavailable.
solver=factor_once(L);
model=struct('m',m,'n',n,'na',na,'dt',dt,'nSteps',nStep,'scheme',scheme, ...
    'nFree',nFree,'A',A,'M',M,'H',H,'W',W,'qX',qX,'qY',qY, ...
    'solve',solver,'R',R,'Ktable',[],'truthFlux',[], ...
    'sourceOrder',sourceOrder);
if buildKernel
    estimatedKernelMB=8*double(nFree)*double(na)*double(nStep+1)/1024^2;
    fprintf('  Precomputing full-order %s flux transfer (about %.1f MB) ...\n',scheme,estimatedKernelMB);
    % K_n^T = dt sum_{j=0}^{n-1} L^{-T}(Q^T)^j H^T; exact recurrence
    % G_0=H^T; tmp_n=L^{-1}G_n; K_{n+1}=K_n+dt*tmp_n;
    % G_{n+1}=R*tmp_n. This avoids the INCORRECT K_{n+1}=Q^T K_n+B.
    tbl=zeros(nFree,na,nStep+1,'double');
    accum=zeros(nFree,na);G=full(H');
    for k=1:nStep
        tmp=solver(G);
        accum=accum+dt*tmp;
        tbl(:,:,k+1)=accum;
        G=R*tmp;
    end
    model.Ktable=tbl;
    clear tbl; % release the temporary workspace reference after transfer
end
fprintf('  FEM ready: DOF=%d, quadrature=%d, %.2fs\n',nFree,nQuad,toc);
end

function solve=factor_once(L)
% Factorization is reused in the 400-step forward/adjoint calculation.
try
    D=decomposition(L,'chol');
    solve=@(B) D\B;
catch
    [LL,UU,p,q]=lu(L,'vector');
    solve=@(B) solve_lu(LL,UU,p,q,B);
end
end

function X=solve_lu(LL,UU,p,q,B)
X=zeros(size(B));
X(q,:)=UU\(LL\B(p,:));
end

function f=create_load(model,kind,xi,b)
mask=reviewer2_shapes('mask',kind,xi,model.qX,model.qY);
f=b*(model.W*mask);
end


function model=build_truth(model,kind,xi,b)
assert(reviewer2_shapes('valid',kind,xi,.001),'Invalid true geometry');
F=create_load(model,kind,xi,b);
U=zeros(model.nFree,1);
model.truthFlux=zeros(model.na,model.nSteps+1);
for k=1:model.nSteps
    U=model.solve(model.R*U+model.dt*F);
    model.truthFlux(:,k+1)=model.H*U;
end
end

function row=angle_interpolate_matrix(vals,theta)
% vals: [na x nVar]; returns [1 x nVar] at continuous angle.
na=size(vals,1);dx=2*pi/na;
x=mod(theta,2*pi)/dx;
el=floor(x/2);z=x-(2*el+1);
w=[.5*z*(z-1),1-z*z,.5*z*(z+1)];
ids=mod(2*el+(0:2),na)+1;
row=w*double(vals(ids,:));
end

function vals=evaluate_truth(model,angles,times)
assert(~isempty(model.truthFlux),'Truth flux was not computed');
times=times(:);angles=angles(:);N=numel(times);
assert(numel(angles)==N,'Angles and times must have the same length');
vals=zeros(N,1);
for j=1:N
    tt=max(0,min(model.nSteps,times(j)/model.dt));
    lo=floor(tt);hi=min(lo+1,model.nSteps);w=tt-lo;
    nflux=(1-w)*model.truthFlux(:,lo+1)+w*model.truthFlux(:,hi+1);
    vals(j)=angle_interpolate_matrix(nflux,angles(j));
end
end

function G=create_design(model,angles,times)
assert(~isempty(model.Ktable),'Backward Euler transfer kernel not available');
times=times(:);angles=angles(:);N=numel(times);
assert(numel(angles)==N,'Angles and times must have the same length');
G=zeros(N,model.nFree);
na=model.na;dx=2*pi/na;
for j=1:N
    tt=times(j)/model.dt;
    if tt<0 || tt>model.nSteps+1e-10,error('Observation time is out of range');end
    lo=min(floor(tt),model.nSteps);hi=min(lo+1,model.nSteps);at=tt-lo;
    a=mod(angles(j),2*pi)/dx;el=floor(a/2);z=a-(2*el+1);
    ww=[.5*z*(z-1),1-z*z,.5*z*(z+1)];
    ia=mod(2*el+(0:2),na)+1;
    Klo=double(model.Ktable(:,ia,lo+1));
    Khi=double(model.Ktable(:,ia,hi+1));
    G(j,:)=(((1-at)*Klo+at*Khi)*ww').';
end
end


function [gx,gy,gw]=source_rectangle_gauss(q)
% Gauss-Legendre nodes for source load only; q=6 exactly original.
if q==6
    [gx,gy,gw]=gaussian_nodes_weights_rect_ref(6);
    return;
end
assert(q>=2 && q==round(q));
kk=(1:q-1)'; off=kk./sqrt(4*kk.^2-1);
J=diag(off,1)+diag(off,-1);
[V,D]=eig(J);[ab,ix]=sort(diag(D));
w=2*(V(1,ix)').^2;
% gx=radial coordinate; gy=angular coordinate; same ordering as original.
[GX,GY]=meshgrid(ab,ab);
gx=GX(:);gy=GY(:);
[Wx,Wy]=meshgrid(w,w);gw=(Wx(:).*Wy(:))';
end
