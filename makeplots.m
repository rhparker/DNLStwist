%% plot 

figure('DefaultAxesFontSize',24);
set(gca,'fontname','times');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

hold on;
lS = {'-','--',':','-.'};
lw = 3;

load dark15
plot(phi,a4,'-','LineWidth',lw);

load dark20
plot(phi,a4,'--','LineWidth',lw);

load dark25
plot(phi,a4,':','LineWidth',lw);

load dark30
plot(phi,a4,'-.','LineWidth',lw);

legendCell = {'$k = 0.15$','$k = 0.20$','$k = 0.25$','$k = 0.30$'};
legend(legendCell,'Interpreter','latex','location','southeast');
xlabel('$\phi$','Interpreter','latex');
ylabel('$a_4$','Interpreter','latex');

set(gca,'xtick',[0:pi/6:pi/3]) % where to set the tick marks
set(gca,'xticklabels',{'$0$','$\pi/6$','$\pi/3$'});

%% k0 vs N for omega = 1

N = [6 8 10 12 14 16 18 20 ...
    22 24 26 28 30 ...
    32 34 36 38 40 ...
    42 44 46 48 50];
k0 = [5.77350E-01 5.41196E-01 5.25731E-01 5.17638E-01 5.12858E-01 5.09796E-01 5.07713E-01 5.06233E-01 ...
    5.05142E-01 5.04314E-01 5.03672E-01 5.03164E-01 5.02754E-01 ...
    5.02419E-01 5.02142E-01 5.01910E-01 5.01714E-01 5.01546E-01 ...
    5.01402E-01 5.01277E-01 5.01168E-01 5.01073E-01 5.00989E-01];

figure('DefaultAxesFontSize',30);
set(gca,'fontname','times');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
ms = 40;

plot(N,k0,'.','MarkerSize',ms);

xlabel('$N$','Interpreter','latex');
ylabel('$k_0$','Interpreter','latex');

%% k0 vs omega for N = 50

figure('DefaultAxesFontSize',30);
set(gca,'fontname','times');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

ms = 40;
lw = 2;

d = -1;
omega = [0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2]
k0 = [0.12524714663, 0.25049429326, 0.37574143989, 0.50098858651, 0.62623573316, ...
    0.7514828798, 0.87673002633, 1.0019764545]
p = polyfit(omega, k0, 1);
hold on;

plot(omega, k0, '.','MarkerSize',ms);
plot(omega, polyval(p,omega),'LineWidth',lw);

xlabel('$\omega$','Interpreter','latex');
ylabel('$k_0$','Interpreter','latex');