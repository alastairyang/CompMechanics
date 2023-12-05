function [U, x] = finiteElementSolve(xl, hs, int_method)
%FINITEELEMENTSOLVE use finite element (1d hat function element) to solve
%the poisson equation in 1D (i.e., -u_xx = f)

%   Input:
%       xl [double]: location of left boundary
%       hs [double array]: Nx1 array of element size, thus location of
%                   right boundary = xl + sum(hs)
%       int_method [string]: numerical integration method. 
%                   'trapz': trapezoidal method (linear interpolation)
%                   'gl': gauss-legendre quadrature (polynomial interp)
%
%   Output:
%       U [double array]: Nx1, solution
%       x [double array]: Nx1, coordinate for plotting

    % basis numbering system
    N = length(hs); % number of elements
    n = (N+1)-1; % number of nodes; minus 1 due to periodic BC 
    K = zeros(n, n);
    hs = [hs; hs(end-1)];

    % stiffness matrix [K]
    for ii = 1:n % iterate over node i
        elmt_i = ii; % 

        nodes_elmt_i = [elmt_i-1,elmt_i, elmt_i+1];
        for jj = 1:n
            if ii == jj % diagonal
                if ii == 1
                    K(ii,jj) = 1./hs(elmt_i) + 1./hs(end);
                elseif ii == n
                    K(ii,jj) = 1./hs(elmt_i) + 1./hs(1);
                else
                    K(ii,jj) = 1./hs(elmt_i) + 1./hs(elmt_i-1);
                end
                % off diagonal
            elseif ismember(ii,nodes_elmt_i) && ismember(jj,nodes_elmt_i)
                K(ii,jj) = -1./hs(elmt_i);
            end

        end
    end

    % forcing vector [f]
    twoL = sum(hs(1:end-1));
    c = -sqrt(pi)./(twoL);
    x = xl + cumsum([0; hs(1:end-1)]);
    x = x(1:end-1); % truncate out the end node bc periodic BC
    % define forcing function
    func = @(x,h,xl,phi) h*phi.*(exp(-(x*h+xl).^2) + c);

    f = zeros(n,1);
    for ii = 1:n
        xe = 0:0.1:1; % local equally spaced coordinate
        % two types of basis in hat element
        phi1 = 1 - xe;
        phi2 = xe; 

        % find the corresponding element length and left node global
        % coordinate value for numerical integration
        if ii == 1 % first node
            xl_phi1 = x(1); % first element
            xl_phi2 = x(end); % last element
            h_phi1 = hs(1); % spacing
            h_phi2 = hs(end);
        else % other interior nodes
            xl_phi1 = x(ii);
            xl_phi2 = x(ii-1);
            h_phi1 = hs(ii); % spacing
            h_phi2 = hs(ii-1);
        end
        % evaluate integral
        switch int_method
            case 'trapz' % trapozoidal
                rhs_phi1 = trapz(xe, func(xe,h_phi1,xl_phi1,phi1));
                rhs_phi2 = trapz(xe, func(xe,h_phi2,xl_phi2,phi2));
            case 'gl' % gauss-legendre
                orderTrunc = 5; % order of truncation for legendre polyn.
                [xx, ww] = lgwt(orderTrunc, xe(1), xe(end));
                rhs_phi1 = sum(func(xx,h_phi1,xl_phi1,xx).*ww);
                rhs_phi2 = sum(func(xx,h_phi2,xl_phi2,fliplr(xx)).*ww);
            otherwise
                error('Unknown numerical integration scheme!')
        end
        % sum contributions from two basis
        f(ii) = rhs_phi1 + rhs_phi2;
    end

    % solve
    U = transpose(K\f);
    x = x';
end
