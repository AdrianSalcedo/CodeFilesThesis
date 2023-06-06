% Function to display probabilities in 3D from t0 to t10
% Jinyan Liu
% Inria Saclay and CMAP Ecole Polytechique
% 2016-2017

clear all;
close all;
y = [-1
-0.9
-0.8
-0.7
-0.6
-0.5
-0.4
-0.3
-0.2
-0.1
0
0.1
0.2
0.3
0.4
0.5
0.6
0.7
0.8
0.9
1];

for k = 0:140
    s1 = 'processLaw/stateDistribution/distributionProba_mode0.t';
    s2 = int2str(k);
    s3 = 'proba at t';
    s4 = 'proba_t';
    filename = strcat(s1,s2);
    figname = strcat(s3,s2);
    savename = strcat(s4, s2);
    x = importdata(filename, '\n');
    disp(x);
    v = reshape(x(:,2), 21, 21);
    figure('name',figname);
    p=pcolor(y,y,v);
    set(p,'edgecolor','none');
    title(figname);
    xlabel('State 1');
    ylabel('State 2');
    c = colorbar;
    colormap (jet(1024));
    ylabel(c,'Density of probability');
    caxis([0,0.013]);  
    % save figures
    %saveas(p,savename,'jpg')
end   
