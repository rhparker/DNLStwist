load evenhole50data

load assym6data

%%  solution

figure('DefaultAxesFontSize',24);
set(gca,'fontname','times');
hold on
lS = {'-','--',':','-.'};
NPlot=4;
for index = 1:NPlot
    plot(t,abs(u(index,:)),'Linewidth',2, 'LineStyle', lS{index} );
end

legendCell = strcat('n=', string(num2cell(1:NPlot)) );legend(legendCell);
xlabel('$z$','Interpreter','latex');
ylabel('$|c_n|$','Interpreter','latex');

%% solution, amplitudes, and phases

figure('DefaultAxesFontSize',24,'Position', [0 0 700 600]);

set(gca,'fontname','times');

ax1 = subplot(2,2,[1 3]);
hold on
lS = {'-','--',':','-.'};
NPlot=4;
for index = 1:NPlot
    plot(t,abs(u(index,:)),'Linewidth',3, 'LineStyle', lS{mod(index,length(lS))+1} );
end
legendCell = strcat('n=', string(num2cell(1:NPlot)) );
legend(legendCell,'Interpreter','latex','location','northeast');
xlabel('$z$','Interpreter','latex');
ylabel('$|c_n|$','Interpreter','latex');

ax2=subplot(2,2,2);
hold on;
plot(1:N,amps,'.b','MarkerSize',40);
% plot(1:N,amps,'-k');
xlabel('$n$','Interpreter','latex');
ylabel('$a_n$','Interpreter','latex');
set(gca,'XTick',1:N);
axis([1,N,-1,1]);
% axis(ax2,'tight');

ax3=subplot(2,2,4);
hold on;
plot(1:N,p,'.r','MarkerSize',40);
% plot(1:N,amps,'-k');
xlabel('$n$','Interpreter','latex');
ylabel('$\theta_n$','Interpreter','latex');
set(gca,'XTick',1:N);
set(gca,'YTick',[-pi/6 0 pi/6]);
axis([1,N,-pi/6 pi/6]);
% yticks([-pi/6 0 pi/6]);
yticklabels({'$-\pi/6$','$0$','$\pi/6$'});
% axis(ax3,'tight');

%% spectrum

figure('DefaultAxesFontSize',24);
set(gca,'fontname','times');
plot(l, '.', 'MarkerSize',40);
axis([-1,1,-2,2]);
xlabel('Re $\lambda$','Interpreter','latex');
ylabel('Im $\lambda$','Interpreter','latex');