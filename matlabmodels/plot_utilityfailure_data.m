close all

load('data/utilityfailurescan','plabel','p','clustersize_distribution','outtage_distribution','iii','jjj','p_join','p_outtage','p_recover','N');

x = p_outtage;
for i = 1:size(clustersize_distribution,2)
    y(i) = find(clustersize_distribution(1,i,:)>0,1,'last')/N^2;
end
plot(x,y)

xlabel('$p_{out}$','Interpreter','latex');
ylabel('cluster size','Interpreter','latex');
