function [gx, gy, gw] = gaussian_nodes_weights_rect_ref(n_gnodes)
%% USAGE: generate Gaussian nodes and weights on the reference rectangular [-1,1]*[-1,1]
% 
% INPUTS:
% n_nodes --- the number of gauss nodes
%             possible values = 1,2,3,4,5,6
%
% OUTPUTS:
% gx --- x coordinates of Gaussian nodes 
% gy --- y coordinates of Gaussian nodes
% gw --- Gaussian weights (row vector)

% Last Modified: 20/12/2016 by zhangqian 
%%
[g,gw] = gaussian_nodes_weights_line_ref(n_gnodes);

[gx, gy]=meshgrid(g, g);

gx=reshape(gx,[],1);

gy=reshape(gy,[],1);

[temp1,temp2]=meshgrid(gw,gw);

temp1=reshape(temp1,[],1);

temp2=reshape(temp2,[],1);

gw=(temp1.*temp2)';