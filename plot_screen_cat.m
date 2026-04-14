function plot_screen_cat(N,name,Nmin,Nmax,cats)

% N1 - growth
% N2 - dynamical resistance
% N3 - steady-state resistance

figure;

h = gscatter(N(:,2),N(:,1),cats,'ygrbk','.',50,'off');
set(h(1), 'Color', '#ffa700'); 
%set(gca,'XScale','log','YScale','log')

hold on
plot([Nmin Nmax],[Nmin Nmax],'--k')

axis([Nmin Nmax Nmin Nmax])

xlabel('Steady-state resistance (\mug/ml)','FontSize',20)
ylabel('Dynamic resistance (\mug/ml)','FontSize',20)
%xticks([0.01 0.1 1 5 10 40 100 1000])
%yticks([0.01 0.1 1 5 10 40 100 1000])
% if strcmp(name, 'spec')
% 
%     ax = ancestor(h, 'axes');
%     ax{5}.YTickLabel = {'1000','100'};
%     ax{5}.XTickLabel = {'1000','100'};
% end


set(findall(gcf,'-property','FontSize'),'FontSize',35)

iptsetpref('ImshowBorder','tight');
set(gcf,'Position',[202   47   920   737])
saveas(gcf,[name '.png'])
end