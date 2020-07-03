load data
g=squeeze(data(4,:,:));
figure
for i=1:40
    for j=2:91
        if(g(i,j)<0.7*g(i,j-1))
            g(i,j:end)=2*g(i,j:end);
        end
    end
    subplot(4,10,i)
    v=diff(smooth(log2(g(i,:)),5,'rlowess'));
    plot(v)
    g(i,2:end)=v;
    g(i,1)=0;
    g([14 15 19],2:end)=0;
    hold on
    plot(g(i,2:end),'r')
    plot([1 90],[0 0],'--k')
    plot([19 19],[0 0.5],'--k')
    axis([0 90 0 0.5])
    axis off
end

save g g
