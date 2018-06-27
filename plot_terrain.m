clear, close all;

[p,t] = loadmesh('sfterrain.off')

p = p';

figure, plotmesh(p,t');
rzview('on')
