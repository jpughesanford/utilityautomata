close all

load('data/utilityfailurescan','plabel','p','clusternumber_distribution','clustersize_distribution','outtage_distribution','iii','jjj','p_join','p_outtage','p_plant','p_recover','N');

clustersize_distribution(1,1,:) = 0;
mesh(p_join,p_outtage,(mean(outtage_distribution,3)))
xlabel('$p_{join}$','Interpreter','latex');
ylabel('$p_{outtage}$','Interpreter','latex');
