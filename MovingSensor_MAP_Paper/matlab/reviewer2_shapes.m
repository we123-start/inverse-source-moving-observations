function out=reviewer2_shapes(action,kind,varargin)
% Original paper's source geometries, with explicit admissibility checks.
% action='truth'|'amplitude'|'physical'|'valid'|'mask'|'boundary'.
kind=lower(char(kind));action=lower(char(action));
switch action
    case 'truth'
        switch kind
            case 'circle', out=[.7;pi/2;.2];
            case 'kite', out=[.4;pi/3;.2];
            case 'fourleaf', out=[.4;pi/2;.7];
            case 'peanut', out=[1;0;0;0;.3];
            otherwise, error('Unknown geometry %s',kind);
        end
    case 'amplitude'
        if strcmp(kind,'peanut'),out=10;else,out=50;end
    case 'physical'
        w=varargin{1}(:);
        prior='H1'; if numel(varargin)>=2 && ~isempty(varargin{2}),prior=char(varargin{2});end
        if strcmp(kind,'peanut')
            a=.5; if strcmpi(prior,'H2'),a=.25;end
            out=[w(1);w(2);w(3);a*w(4);a*w(5)];
        else
            out=[atan(w(1))/pi+.5; 2*atan(w(2))+pi; atan(w(3))/pi+.5];
        end
    case 'valid'
        xi=varargin{1}(:);
        margin=varargin{2};
        out=shape_valid(kind,xi,margin);
    case 'mask'
        xi=varargin{1}(:);x=varargin{2};y=varargin{3};
        out=inside_shape(kind,xi,x,y);
    case 'boundary'
        xi=varargin{1}(:);ang=varargin{2};
        [xx,yy]=boundary_shape(kind,xi,ang);out=[xx(:)';yy(:)'];
    otherwise
        error('Unknown shape action: %s',action);
end
end

function ok=shape_valid(kind,xi,margin)
if any(~isfinite(xi)),ok=false;return;end
if strcmp(kind,'peanut')
    if numel(xi)~=5,ok=false;return;end
    a=linspace(0,2*pi,1001);
    r=xi(1)/2+xi(2)*cos(a)+xi(3)*sin(a)+xi(4)*cos(2*a)+xi(5)*sin(2*a);
    % Certified grid bound: |q'| <= sum_i i*hypot(a_i,b_i), so
    % extrema between grid points deviate by at most Lipschitz*h/2.
    lip=hypot(xi(2),xi(3))+2*hypot(xi(4),xi(5));
    safety=lip*pi/(numel(a)-1);
    ok=min(r)-safety>margin && max(r)+safety<1-margin;
else
    if numel(xi)~=3,ok=false;return;end
    rho=xi(1);phi=xi(2);s=xi(3);
    if ~(rho>0 && rho<1 && phi>0 && phi<2*pi && s>0 && s<1)
        ok=false;return;
    end
    switch kind
        case 'circle'
            ok=rho+s<1-margin;
        case 'fourleaf'
            % r(theta) = s*(0.4 + 0.2*cos(4*theta)); sufficient bound.
            ok=rho+.6*s<1-margin;
        case 'kite'
            t=linspace(0,2*pi,1001);
            [x,y]=boundary_shape('kite',xi,t);
            % Certified bound from ||(x',y')|| <= s*hypot(2.3,1.5).
            lip=s*hypot(2.3,1.5);
            safety=lip*pi/(numel(t)-1);
            ok=max(hypot(x,y))+safety<1-margin;
        otherwise, error('Unknown geometry %s',kind);
    end
end
end

function mask=inside_shape(kind,xi,x,y)
switch kind
    case 'circle'
        xc=xi(1)*cos(xi(2));yc=xi(1)*sin(xi(2));
        mask=((x-xc).^2+(y-yc).^2 <= xi(3)^2);
    case 'fourleaf'
        xc=xi(1)*cos(xi(2));yc=xi(1)*sin(xi(2));
        dx=x-xc;dy=y-yc;a=atan2(dy,dx);
        mask=(hypot(dx,dy) <= xi(3)*(.4+.2*cos(4*a)));
    case 'peanut'
        a=atan2(y,x);
        q=xi(1)/2+xi(2)*cos(a)+xi(3)*sin(a)+xi(4)*cos(2*a)+xi(5)*sin(2*a);
        mask=(hypot(x,y)<=q);
    case 'kite'
        % Analytic inside-test avoids a costly 241-edge inpolygon for every parameter estimate.  In local coordinates v=y/(1.5*s), the two
        % boundary branches are u=-1.3*v^2 +/- sqrt(1-v^2).
        % This follows directly from y=1.5*s*sin(theta),
        % x=s*(cos(theta)+0.65*cos(2theta)-0.65).
        xc=xi(1)*cos(xi(2));yc=xi(1)*sin(xi(2));s=xi(3);
        v=(y-yc)/(1.5*s);u=(x-xc)/s;
        rr=sqrt(max(0,1-v.^2));
        mask=abs(v)<=1 & u>=-1.3*v.^2-rr & u<=-1.3*v.^2+rr;
    otherwise, error('Unknown geometry %s',kind);
end
mask=double(mask(:));
end

function [x,y]=boundary_shape(kind,xi,a)
a=a(:)';
switch kind
    case 'circle'
        r=xi(3)+0*a;
        x=xi(1)*cos(xi(2))+r.*cos(a);
        y=xi(1)*sin(xi(2))+r.*sin(a);
    case 'fourleaf'
        r=xi(3)*(.4+.2*cos(4*a));
        x=xi(1)*cos(xi(2))+r.*cos(a);
        y=xi(1)*sin(xi(2))+r.*sin(a);
    case 'kite'
        % Entire bracket multiplied by scaling factor xi(3), as in source code.
        x=xi(1)*cos(xi(2))+xi(3)*(cos(a)+.65*cos(2*a)-.65);
        y=xi(1)*sin(xi(2))+1.5*xi(3)*sin(a);
    case 'peanut'
        r=xi(1)/2+xi(2)*cos(a)+xi(3)*sin(a)+xi(4)*cos(2*a)+xi(5)*sin(2*a);
        x=r.*cos(a);y=r.*sin(a);
    otherwise, error('Unknown geometry %s',kind);
end
end

