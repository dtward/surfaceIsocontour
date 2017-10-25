% test
clear all;
close all;

% your surface is stored with faces f and vertices v
addpath /cis/home/dtward/Functions/byu
[f,v] = readbyu_local('/cis/home/dtward/Documents/Mostofsky/atlas_2693/atlas_2693_01.byu');

% here's some signal, suppose this is p-values at each vertex
data = ( cos(2*pi*v(:,1)/40) + cos(2*pi*v(:,2)/50) + cos(2*pi*v(:,3)/60) )*3;
pvalues = 1-chi2cdf(abs(data),1);  % some fake stats for this example
h = patch('faces',f,'vertices',v,'facevertexcdata',data,'facecolor','interp','facelighting','gouraud','edgecolor','none');
axis image
view(16,18)
colorbar
drawnow;

% do the isocontour
addpath /cis/home/dtward/Documents/sandbox/surfaceIsocontour
threshold = 0.05;
[p,e] = surfaceIsocontour(f,v,pvalues,threshold);
% this returns vertices p, and edges (aka faces with two elements) e
patch('faces',e,'vertices',p,'facecolor','none','edgecolor','k','linewidth',2);
drawnow;