function [x, omega] = bendingBeamFEMFunc(n)
%BENDINGBEAMFEMFUNC Finite element modeling of cantilever beam deflection

    % constants
    F = 40; % newton 
    q0 = 100; % newton per meter
    L = 0.1; % meter, 100 mm in doc
    E = 210e9; % Young modulus N/m^2
    Iy = 1/12*0.01^4; % area moment

    % build the linear system
    f = zeros(n, 1); % Forcing vector
    K = zeros(n, n); % Stiffness matrix
    % stiffness matrix
    for ii = 2:n
        for jj = 2:n
            if ii+jj ~= 3
                K(ii,jj) = ii*jj*(ii-1)*(jj-1)*(2*L)^(ii+jj-3)/(ii+jj-3);
            else
                K(ii,jj) = 0;
            end
        end
    end
    % forcing vector
    for jj = 2:n
        f(jj) = (1/(E*Iy))*( q0*L^(jj+1)/(jj+1) -...
                             (q0*L^(jj+2)*(2^(jj+2)-1))/(L*(jj+2)) + ...
                             2*q0*L^(jj+1)*(2^(jj+1)-1)/(jj+1) + ...
                             (2*L)^jj*F);
    end
    % truncate the first order terms which are the first row and column
    K = K(2:n,2:n);
    f = f(2:n);
    % solve for coefficient a
    a = K\f;
   
    % construct solution
    x = 0:(2*L)/100:2*L;
    omega = 0;
    for kk = 1:n-1
        p = kk+1;
        omega = omega + a(kk)*x.^p;
    end
    
end

