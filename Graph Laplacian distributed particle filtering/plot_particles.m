function plot_particles(n, t, p, w)
%
% Description:
% 'plot the particle distribution .'
%
% INPUTS:
%  p - locations of nodes (or to say, particles) 
%  w - weights(probability)
%  t - time span
%  n - number of sensors
%
% OUTPUTS:
%  figure
%


w = ceil(1+(w-min(w))/(max(w)-min(w)+0.1)*100);
N = size(w,1);


jet_color = colormap(jet(max(w)));
for i = 1:N
    selected_color = jet_color(w(i),:);
    plot(p(1,i), p(2,i), 'o','color','k','MarkerEdgeColor',selected_color,...
     'MarkerFaceColor',selected_color,'MarkerSize',4,'LineWidth',0.2);  % plot weights by varying colors
    hold on;
end
c.Label.String = 'weight';
title(['Particle Distribution,(t=',num2str(t),',N=',num2str(n),')']);


end