function [c,ceq]=map_constraints(w,kind)
xi=reviewer2_shapes('physical',kind,w,'H1');
if strcmp(kind,'circle')
    c=xi(1)+xi(3)-1+1e-6;
elseif strcmp(kind,'fourleaf')
    c=xi(1)+.6*xi(3)-1+1e-6;
elseif strcmp(kind,'kite')
    a=linspace(0,2*pi,1001);B=reviewer2_shapes('boundary',kind,xi,a);
    safety=xi(3)*hypot(2.3,1.5)*pi/1000;
    c=max(hypot(B(1,:),B(2,:)))+safety-1+1e-6;
else
    a=linspace(0,2*pi,1001);
    q=xi(1)/2+xi(2)*cos(a)+xi(3)*sin(a)+xi(4)*cos(2*a)+xi(5)*sin(2*a);
    % Parameter-space Lipschitz safety, consistent with shape_valid.
    lip=hypot(xi(2),xi(3))+2*hypot(xi(4),xi(5));
    safety=lip*pi/(numel(a)-1);
    c=[1e-6+ safety-min(q); max(q)+safety-1+1e-6];
end
ceq=[];c=c(:);
end
