k = 0.25;
d = -1;

t = linspace(0,4*pi,1000);

% complex version

% % hole N=7
% % mags = [ 8.24135E-01  -2.33386E-01   5.85478E-02 0 -5.85478E-02 2.33386E-01 -8.24135E-01]';
% % double hole N=6
% % mags = [ 0 9.65925E-01  -2.58822E-01 0 2.58822E-01 -9.65925E-01 ]';

load mags;
mags = mags';

% even points, two opposite holes
% mags = [mags -flip(mags(2:end-1))];

% odd points, hole in center
% mags = [mags -flip(mags(1:end-1))];

% generate phases evenly spaces
% N = length(mags);
% phi = 2*pi/N;
% phi = 0; 
% nn = 1:N;
% u0 = exp( 1i*phi*nn');
% u0 = u0.*mags';

% % get phases from AUTO
% N = (length(mags)+1)/2;
% p = [0 ; mags(N+1:end)];
% u0 = ( mags(1:N).*exp( 1i*p ) );
% phi = 0.25;
% % g = 0;
% % phi = pi/N; 

w = 1;

% even hole from AUTO, phi = pi/N
N = length(mags)*2;
phi = pi/N;
nn = [1:length(mags) - 1]';
p = [0 ; nn ; 0 ; -flip(nn) ]*phi;
u0 = [mags ; 0 ; flip(mags(2:end)) ].*exp(1i*p);

% % odd hole from AUTO, phi = pi/N
% N = length(mags)*2 + 1;
% phi = pi/N;
% nn = [1:length(mags)]';
% p = [0 ; nn*phi-pi/2 ; -flip(nn*phi-pi/2) ];
% u0 = [0 ; mags ; flip(mags) ].*exp(1i*p);

% perturbation 
% u0(1) = u0(1) * 1.05;

% k = 0.25;

% solve on interval  with IC (regular)
u  = rk4( @(s,u) twist(s,u,k,phi,d), u0, t);

% solve on interval with IC (g-version)
% g = 0.1;
% u0(4) = u0(4)+.01;
% u0 = [1 0 0 0 0 0]';

% u  = rk4( @(s,u) twist_g(s,u,k,phi,d,g), u0, t);

% % solve on interval with IC (k-version)
% k = 0.2*ones(size(u0));
% k(1) = 0.4;
% phi = 0.25;
% % u0(4) = u0(4)+.01;
% u  = rk4( @(s,u) twist_k(s,u,k,phi,d), u0, t);
% w = 1;

w = 1;
u1 = u0.*exp(1i*w*t);

J = twistJ(real(u0),imag(u0),k,phi,d,w);
[V,l] = eig(J);
l = diag(l);

phases = angle(u0);
amps = u0;
for index = 1:length(u0)
    if (real(u0(index)) < 0)
        amps(index) = -abs( u0(index) );
        phases(index) = angle( u0(index) ) + pi;
    else
        amps(index) = abs( u0(index) );
        phases(index) = angle( u0(index) );
    end
end

% AClimit = 0*u0;
% AClimit(1) = 1;
% J0 = twistJ(AClimit,0*AClimit,0,phi,d,w);
% [V0,l0] = eig(J0);
% l0 = diag(l0);

figure('DefaultAxesFontSize',24);
set(gca,'fontname','times');
plot(t,abs(u),'Linewidth',3 );
% plot(t,abs(u).^2,'Linewidth',3 );
legendCell = string(num2cell(1:N));
legend(legendCell);
xlabel('$z$','Interpreter','latex');
ylabel('$|c_n|$','Interpreter','latex');

figure('DefaultAxesFontSize',20);
set(gca,'fontname','times');
% spectrum plot
plot(l, '.', 'MarkerSize',30);
axis([-1e-12,1e-12,-2,2]);
xlabel('Re $\lambda$','Interpreter','latex');
ylabel('Im $\lambda$','Interpreter','latex');


%%

% make plots

figure('DefaultAxesFontSize',24,'Position', [0 0 1600 600]);
set(gca,'fontname','times');
subplot(1,2,1);
plot(t,real(u),'Linewidth',3 );
legendCell = string(num2cell(1:N));
legend(legendCell,'Interpreter','latex');
xlabel('$z$','Interpreter','latex');
ylabel('Re $c_n$','Interpreter','latex');

subplot(1,2,2);
hold on;
plot(1:N,amps,'.','MarkerSize',30);
plot(1:N,amps,'-k');
xlabel('$n$','Interpreter','latex');
ylabel('$a_n$','Interpreter','latex');
% 
% subplot(1,2,2);
% hold on;
% plot(1:N,abs(amps),'.','MarkerSize',30);
% plot(1:N,abs(amps),'-k');
% xlabel('$n$','Interpreter','latex');
% ylabel('$|c_n|$','Interpreter','latex');


%% 
d = -1;
omega = [0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2]
k0 = [0.12524714663, 0.25049429326, 0.37574143989, 0.50098858651, 0.62623573316, ...
    0.7514828798, 0.87673002633, 1.0019764545]
p = polyfit(omega, k0, 1);
figure('DefaultAxesFontSize',20);
hold on;
set(gca,'fontname','times');
plot(omega, k0, '.','MarkerSize',30);
plot(omega, polyval(p,omega),'LineWidth',2);
xlabel('$\omega$','Interpreter','latex');
ylabel('$k_0$','Interpreter','latex');


%%

% plot bifurcation diagram
load bd;
figure('DefaultAxesFontSize',20);
set(gca,'fontname','times');
plot(par,l2norm,'LineWidth',3);
xlabel('$k$','Interpreter','latex');
ylabel('$l^2$ norm','Interpreter','latex');


%% try shooting method

t = linspace(0,2*pi,1000);
phi = pi/8;
g = 0.1;
k = 0.25;
options = optimoptions('fsolve', 'Algorithm','levenberg-marquardt','Display','iter');
a = [1 0 0 0 0 0]';
a1 = fsolve(@(y) BC_g(t, y, k, phi, d, g), a, options);
u = rk4( @(s,u) twist_g(s,u,k,phi ,d, g), a1, t);
%%
figure;
subplot(1,2,1);
plot(t,real(u),'Linewidth',3);
subplot(1,2,2);
plot(t,abs(u),'Linewidth',3);

%%

% phases

% make everything have positive real part (we allow coeffs to be negative)

a1 = u0;
m1 = u0;
phases = angle(a1);

for index = 1:length(a1)
    if (real(a1(index)) < 0)
        m1(index) = -abs( a1(index) );
        phases(index) = angle( a1(index) ) + pi;
    else
        m1(index) = abs( a1(index) );
        phases(index) = angle( a1(index) );
    end
end



%% real version

% nn = 0:5;
% a0 = cos(phi*nn');
% b0 = sin(phi*nn');
% u = rk4( @(s,u) twistreal(s,u,k,phi,d), [a0 ; b0], t);
% plot(t, u(1:6,:) );

%% subroutines
 
function u = BC(t, u0, k, phi, d)
    u1 = rk4( @(s,u) twist(s,u,k,phi,d), u0, t);
    u = [u1(:,end) - u0];
end

function u = BC_g(t, u0, k, phi, d, g)
    u1 = rk4( @(s,u) twist_g(s,u,k,phi,d,g), u0, t);
    u = [u1(:,end) - u0];
end

function u = BC_t(t, u0, k, phi, d)
    u1 = rk4( @(s,u) twist_t(s,u,k,phi,d), u0, t);
    u = [u1(:,end) - u0];
end

function dudt = twist(t,u,k,phi,d)
    N = length(u);
    K =  exp(-1i*phi)*diag( ones(1,N-1), 1  ) ...
        + exp(1i*phi)*diag( ones(1,N-1), -1 );
    K(1,N) = exp(1i*phi);
    K(N,1) = exp(-1i*phi);
    Nc = diag( abs(u).^2 );
    dudt = -1i*(k*K*u + d*Nc*u);
end

function dudt = twist_w(t,u,k,phi,d,w)
    N = length(u);
    K =  exp(-1i*phi)*diag( ones(1,N-1), 1  ) ...
        + exp(1i*phi)*diag( ones(1,N-1), -1 );
    K(1,N) = exp(1i*phi);
    K(N,1) = exp(-1i*phi);
    Nc = diag( abs(u).^2 );
    dudt = -1i*(k*K*u + w*eye(N)*u + d*Nc*u);
end

% PT symmetric version with gain/loss
function dudt = twist_g(t,u,k,phi,d,g)
    N = length(u);
    K =  exp(-1i*phi)*diag( ones(1,N-1), 1  ) ...
        + exp(1i*phi)*diag( ones(1,N-1), -1 );
    K(1,N) = exp(1i*phi);
    K(N,1) = exp(-1i*phi);
    Nc = diag( abs(u).^2 );
    GL = diag( (-1).^(2:N+1) );
    dudt = -1i*(k*K*u + d*Nc*u + 1i*g*GL*u);
end

% PT symmetric version with gain/loss and omega
function dudt = twist_g_w(t,u,k,phi,d,g,w)
    N = length(u);
    K =  exp(-1i*phi)*diag( ones(1,N-1), 1  ) ...
        + exp(1i*phi)*diag( ones(1,N-1), -1 );
    K(1,N) = exp(1i*phi);
    K(N,1) = exp(-1i*phi);
    Nc = diag( abs(u).^2 );
    GL = diag( (-1).^(2:N+1) );
    dudt = -1i*(k*K + d*Nc + w*eye(N) + 1i*g*GL)*u;
end

% in this case k is a vector, so can have different couplings
function dudt = twist_k(t,u,k,phi,d)
    N = length(u);
    K =  exp(-1i*phi)*diag( k(1:end-1), 1  ) ...
        + exp(1i*phi)*diag( k(1:end-1), -1 );
    K(1,N) = exp(1i*phi)*k(end);
    K(N,1) = exp(-1i*phi)*k(end);
    Nc = diag( abs(u).^2 );
    dudt = -1i*(K*u + d*Nc*u);
end

function dudt = twist_t(t,u,k,phi,d)
    N = length(u);
    K =  exp(-1i*phi)*diag( ones(1,N-1), 1  ) ...
        + exp(1i*phi)*diag( ones(1,N-1), -1 );
    K(1,N) = exp(1i*phi*t);
    K(N,1) = exp(-1i*phi*t);
    Nc = diag( abs(u).^2 );
    dudt = -1i*(k*K*u + d*Nc*u);
end

function dudt = twistreal(t,u,k,phi,d)
    N = length(u)/2;
    a = u(1:N);
    b = u(N+1:end);
    dadt = [];
    dbdt = [];
    for r = 1:N
        % indices for periodic BCs
        prev = perindex(r-1,N);
        next = perindex(r+1,N);
        dadt = [dadt ;  k*( -sin(phi)*a(next) +cos(phi)*b(next) ...
                            +sin(phi)*a(prev) +cos(phi)*b(prev) ) ...
                            + d*(a(r)^2 + b(r)^2)*b(r) ];
        dbdt = [dbdt ; -k*(  cos(phi)*a(next) +sin(phi)*b(next) ...
                            +cos(phi)*a(prev) -sin(phi)*b(prev) ) ...
                            - d*(a(r)^2 + b(r)^2)*a(r) ];                          
    end
    dudt = [dadt ; dbdt];
end

% Jacobian for twisted system
function J = twistJ(a,b,k,phi,d,w)
    N = length(a);
    Id = eye(N);
    Z = zeros(N,N);
    wblock = w * [ [ Z Id ] ; [ -Id Z ] ];
    NLblock = d * [ [ diag(2*a.*b) diag(a.^2+3*b.^2) ] ; ...
                    [ diag(-(3*a.^2+b.^2)) diag(-2*a.*b) ] ];
    UD1 = diag( ones(1,N-1), 1  ); UD1(N,1) = 1;
    LD1 = diag( ones(1,N-1), -1 ); LD1(1,N) = 1;
    C = cos(phi)*(  UD1 + LD1 );
    S = sin(phi)*( -UD1 + LD1 );
    kblock = k* [ [ S C ] ; [ -C S ] ];
    J = kblock + wblock + NLblock;
end

% index for periodic BCs
function n = perindex(r,N)
    if r == 0
        n = N;
    elseif r == N+1
        n = 1;
    else
        n = r;
    end
end

% Runge-Kutta 4 ODE solver
% t is time grid
function u = rk4(f, u0, t)
    u = u0;
    h = t(2) - t(1);
    for index = 1:(length(t) - 1)
       k1 = h*f( t(index), u(:,end) );
       k2 = h*f( t(index)+h/2, u(:,end)+0.5*k1 );
       k3 = h*f( t(index)+h/2, u(:,end)+0.5*k2 ); 
       k4 = h*f( t(index)+h, u(:,end)+k3 );
       u = [ u  u(:,end)+(k1 + 2*k2 + 2*k3 + k4)/6 ];
    end
end
