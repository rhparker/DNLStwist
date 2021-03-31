% load oddhole7pert
% load evenhole6pert
% load evenhole10k045pertk035
load evenhole10k045pertk055
% load dp12k025pert
% load assym6pert


%%  solution

figure('DefaultAxesFontSize',40);
set(gca,'fontname','times');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on
lS = {'-','--',':','-.'};
NPlot=4;
for index = 1:NPlot
    plot(t,abs(u(index,:)),'Linewidth',3, 'LineStyle', lS{index} );
end

legendCell = strcat('n=', string(num2cell(1:NPlot)) );legend(legendCell);
xlabel('$z$','Interpreter','latex');
ylabel('$|c_n|$','Interpreter','latex');

% axis([0 20*pi 0 1]);
axis([0 40*pi 0 0.6]);
% set(gca,'XTick', [0 25*pi 50*pi 75*pi 100*pi] );
% xticklabels({'$0$','$25\pi$','$50\pi$','$75\pi$','$100\pi$'});
% set(gca,'XTick', [0 10*pi 20*pi] );
% xticklabels({'$0$','$10\pi$','$20\pi$'});
set(gca,'XTick', [0 20*pi 40*pi] );
xticklabels({'$0$','$20\pi$','$40\pi$'});


