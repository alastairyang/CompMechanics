function [U, x] = spectralSolve(nx, totL)
%SPECTRALSOLVE Solving 1D poisson equation with periodic BC (i.e., -u_xx = f) 
% using spectral method
% Input:
%       nx [double]: number of sampled spacing points
%       totL [double]: total length of the domain, i.e., 2L
% Output:
%       


    if mod(nx,2) ~= 0
        error('Number of wavenubmer has to be an even number!')
    end

    % construct spatial grid 
    x = linspace(-totL/2, totL/2, nx);
    func = @(x) exp(-x.^2) - sqrt(pi)/(totL);
        
    % wavenumber k
    N = length(x);
    ni = -N/2:N/2-1 ;
    kn = ni*2*pi/totL;
    id_zero = find(ni==0);
    Frho = fftshift(fft(func(x)));
    % remove k=0 case (as its corresponding \rho is 0 anyways)
    Frho(id_zero) = [];
    kn(id_zero) = [];
    Fu = Frho./(kn.^2);
    % add zero wavenumber case back
    Fu = [Fu(1:id_zero-1) 0 Fu(id_zero:end)];

    % shift back, and take inverse FFT
    U = ifft(ifftshift(Fu));

end

