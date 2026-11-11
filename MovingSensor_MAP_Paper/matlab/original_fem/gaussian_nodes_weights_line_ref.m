function [g,gw] = gaussian_nodes_weights_line_ref(n_gnodes)
%% USAGE: generate Gaussian nodes and weights on the reference segment [-1,1]
%
% INPUTS:
% n_nodes --- the number of gauss nodes
%             possible values = 1,2,3,4,5,6
%
% OUTPUTS:
% g --- x coordinates of Gaussian nodes 
% gw --- Gaussian weights (row vector)
% Last Modified: 28/12/2016 by Qian Zhang


%%
if n_gnodes == 1
 
    gw = 2;
    g=0;
 
elseif n_gnodes == 2
 
    gw = [1,1];
    g(1) = -0.57735026918962576451;
    g(2) = 0.57735026918962576451;
 
elseif n_gnodes == 3
 
    gw=[5/9 8/9 5/9];
    g(1) = -0.77459666924148337704;
    g(2) = 0;
    g(3) = 0.77459666924148337704;
 
elseif n_gnodes == 4
 
    gw = [0.347854845137454, 0.652145154862546, ...
          0.652145154862546, 0.347854845137454];
    g(1) = -sqrt(3+2*sqrt(6/5))/sqrt(7);
    g(2) = -sqrt(3-2*sqrt(6/5))/sqrt(7);
    g(3) = -g(2);
    g(4) = -g(1);
 
elseif n_gnodes == 5
 
    gw = zeros(1,n_gnodes);
    gw(1) = 0.2369268850561890875;
    gw(2) = 0.478628670499366468;
    gw(3) = 0.56888888888888888889;
    gw(4) = gw(2);
    gw(5) = gw(1);
    
    g(1) = -0.9061798459386639928;
    g(2) = -0.53846931010568309104;
    g(3) = 0;
    g(4) = -g(2);
    g(5) = -g(1);
 
elseif n_gnodes == 6
    
    gw = zeros(1,n_gnodes);
    gw(1) = 0.17132449237917034504;
    gw(2) = 0.36076157304813860757;
    gw(3) = 0.46791393457269104739;
    gw(4) = gw(3);
    gw(5) = gw(2);
    gw(6) = gw(1);
    
    g(1) = -0.93246951420315202781;
    g(2) = -0.66120938646626451366;
    g(3) = -0.23861918608319690863;
    g(4) = 0.23861918608319690863;
    g(5) = -g(2);
    g(6) = -g(1);
    
else
    
    disp('number of Gaussian nodes should be between 1 and 6.')
    
end
 
return;