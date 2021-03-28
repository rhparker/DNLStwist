k = 0.25;
d = -1;
w = 1;

t = linspace(0,100*pi,20000);

load mags6;
mags = mags';
N = 6;

% get phases from AUTO
p = mags(N+1:end) - mags(N+1);
amps = mags(1:N);

% perturbation
amps(4) = amps(4)+0.05;

u0 = amps.*exp( 1i*p );
phi = 0.25;
w = 1;
k = [0.4 0.25 0.25 0.25 0.25 0.25];



% solve on interval
u  = rk4( @(s,u) twist_k(s,u,k,phi,d), u0, t);

J = twistJ_k(real(u0),imag(u0),k,phi,d,w);
[V,l] = eig(J);
l = diag(l);

% solution
figure('DefaultAxesFontSize',24);
set(gca,'fontname','times');
hold on
lS = {'-','--',':','-.'};
NPlot=4;
for index = 1:NPlot
    plot(t,abs(u(index,:)),'Linewidth',2, 'LineStyle', lS{index} );
end
% plot(t,abs(u).^2,'Linewidth',3 );
legendCell = strcat('n=', string(num2cell(1:NPlot)) );legend(legendCell);
xlabel('$z$','Interpreter','latex');
ylabel('$|c_n|$','Interpreter','latex');

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

% Jacobian for twisted system, when k is a vector
function J = twistJ_k(a,b,k,phi,d,w)
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
    kmask = diag(k(1:N-1),1) + diag(k(1:N-1),-1);
    kmask(1,N) = k(N);
    kmask(N,1) = k(N);
    kblock = [ [ kmask.*S kmask.*C ] ; [ -kmask.*C kmask.*S ] ];
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
