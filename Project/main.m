%% Computational Methods in Mechanics: Final project
% Author: Donglai Yang
% Dec 5, 2023
%% global parameter
totL = 20; % total length of domain, 2L
xl = -10;
nx = 256;
nk = 48; % No. wavenumbers
hs = 0.1*ones(200,1);
h = hs(1);

%% Finite difference
maxlogN = 3;
N = 10^maxlogN;
M = 10; % number of samples
[U_fd, x_fd] = finiteDifferenceSolve(xl, totL, N);

figure;
plot(x_fd, U_fd)
xlabel('x'); ylabel('u')
title('Finite difference solution to -u_{xx} = \rho')

% check convergence order
Ns = linspace(10,N,M);

Es = zeros(length(Ns),1);
for ii = 1:length(Ns)
    [U_fd, x_fd] = finiteDifferenceSolve(xl, totL, Ns(ii));
    Es(ii) = computeEnergy(x_fd, U_fd, totL);
end

figure;
Eref = Es(end);
loglog(Ns, abs(Es - Eref));
xlabel('Number of elements','FontSize',14)
ylabel('Residue in energy', 'FontSize',14)

%% Finite element
% numerical integration: gauss-legendre
int_method = 'gl';
[U_fe, x_fe] = finiteElementSolve(xl, hs,int_method);

figure;
plot(x_fe, U_fe)
xlabel('x'); ylabel('u')
title('Finite element solution to -u_{xx} = \rho')

% test convergence
maxlogN = 3;
N = 10^maxlogN;
M = 10; % number of samples
Ns = linspace(10,N,M);
Es = zeros(length(Ns),1);
for ii = 1:length(Ns)
    hs = totL/Ns(ii)*ones(Ns(ii),1);
    [U_fe, x_fe] = finiteElementSolve(xl, hs, int_method);
    Es(ii) = computeEnergy(x_fe, U_fe, totL);
end

figure;
Eref = Es(end);
loglog(Ns, abs(Es - Eref));
xlabel('Number of elements','FontSize',14)
ylabel('Residue in energy', 'FontSize',14)

%% Spectral method
[U_st, x_st] = spectralSolve(nx, totL);

figure;
plot(U_st)
xlabel('x'); ylabel('u')
title('Spectral method solution to -u_{xx} = \rho')

% test convergence
maxlogN = 14;
Ns = 2.^(3:maxlogN);
Es = zeros(length(Ns),1);
for ii = 1:length(Ns)
    [U_st, x_st] = spectralSolve(Ns(ii), totL);
    Es(ii) = computeEnergy(x_st, U_st, totL);
end

figure;
Eref = Es(end);
loglog(Ns, abs(Es - Eref));
xlabel('Number of elements','FontSize',14)
ylabel('Residue in energy', 'FontSize',14)

%% Appendix: functions
function E = computeEnergy(x,u, totL)
    % the energy is \int_{-L}^{L} -(1/2)*(du/dx)^2 + \rho*u dx

    % define forcing function
    rhoFunc = @(x) exp(-x.^2) - sqrt(pi)/(totL);

    dudx = computeDuDx(x,u);
    neg_dudx_sqr = -1*dudx.^2;

    fu = rhoFunc(x).*u;
    
    % use Trapezoidal rule
    E = trapz(x, 0.5*neg_dudx_sqr + fu);
end

function dudx = computeDuDx(x, u)
    % center difference in the middle and 2nd order forward/backward diff
    % on the boundary
    h = x(2)-x(1); % assuming uniform spacing
    dudx_l = (-3*u(1)+4*u(2)-u(3))/(2*h);
    dudx_r = -(u(end-2)-4*u(end-1)+3*u(end))/(2*h);
    dudx_c = (u(3:end) - 2*u(2:end-1) + u(1:end-2)) ./ (x(3:end) - x(1:end-2));
    dudx = [dudx_l, dudx_c, dudx_r];
    
end

