function err = calcEnergyNorm(w_exact, w_approx, x, q)
%CALCENERGYNORM Summary of this function goes here
    F = 40; % newton 
    E = 210e9; % Young modulus N/m^2
    Iy = 1/12*0.01^4; % area moment    

    % calculation
    dx = mean(x(2:end) - x(1:end-1));
    w_diff = w_exact - w_approx;
    w_diff_xx = diff(w_diff,2)/(dx^2);
    % due to 2nd order derivative, we need to truncate original derivative
    % size
    integrant = 0.5*E*Iy*w_diff_xx - q(2:end-1).*w_diff(2:end-1);
    integral = trapz(integrant, x(2:end-1));
    forcing = F*w_diff(end);
    err = abs(integral + forcing);


end

