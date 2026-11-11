function [Aa,Athetha,M]=superstiffmatrixpolar(xy,iiU,hx,hy)
global gx gy gw
Jj=hx*hy;
np=size(xy,2);
nt=size(iiU,2);
% A=sparse(np,np);
Aa=sparse(np,np);%r*v_r*v_r,v_r代表的是对r的求导,v是基函数
Athetha=sparse(np,np);%1/r*v_theta*v_theta, v_theta代表的是对theta求导
M=sparse(np,np);
[Phau,Pcau,Peau]=superbasicu(gx',gy',hx,hy);
for Le=1:nt
    PointU=iiU(:,Le);
    xe=xy(1,PointU);
    % ye=xy(2,PointU);
    x0=sum(xe)/9;
    % y0=sum(ye)/9;
    xx=hx*gx+x0*ones(length(gx),1);
    % yy=hy*gy+y0*ones(length(gy),1);
    bbE=repmat(gw,9,1);
    Ur=(bbE.*Pcau.*xx')*(Pcau')*Jj;
    Utheta=(bbE.*Peau./xx')*(Peau')*Jj;
    mUv=(bbE.*Phau.*xx')*(Phau')*Jj;
    % bbE=repmat(gw,9,1);
    % gradUv=((bbE.*Pcau)*(Pcau')+...
    %     (bbE.*Peau)*(Peau'))*Jj;
    % mUv=(bbE.*Phau)*(Phau')*Jj;
    % A(PointU,PointU)=A(PointU,PointU)+gradUv;
     Aa(PointU,PointU)=Aa(PointU,PointU)+Ur;
    Athetha(PointU,PointU)=Athetha(PointU,PointU)+Utheta;

    M(PointU,PointU)=M(PointU,PointU)+mUv;
end