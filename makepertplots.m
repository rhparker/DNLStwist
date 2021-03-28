% load oddhole7pert
% load evenhole6pert
% load evenhole10k045pertk035
% load evenhole10k045pertk055
% load dp12k025pert
load assym6pert


%%  solution

figure('DefaultAxesFontSize',36);
set(gca,'fontname','times');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on
lS = {'-','--',':','-.'};
NPlot=4;
for index = 1:NPlot
    plot(t,abs(u(index,:)),'Linewidth',2, 'LineStyle', lS{index} );
end

legendCell = strcat('n=', string(num2cell(1:NPlot)) );legend(legendCell);
xlabel('$z$','Interpreter','latex');
ylabel('$|c_n|$','Interpreter','latex');

axis([0 100*pi 0 1]);
% set(gca,'XTick', [0 25*pi 50*pi 75*pi 100*pi] );
% xticklabels({'$0$','$25\pi$','$50\pi$','$75\pi$','$100\pi$'});
set(gca,'XTick', [0 50*pi 100*pi] );
xticklabels({'$0$','$50\pi$','$100\pi$'});


