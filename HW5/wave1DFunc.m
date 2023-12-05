function U = wave1DFunc(L, t, dt)
%WAVE1DFUNC simulate 1D wave equation utt = uxx using explicit finite
%difference method

    c = 1;
    dx = dt/c;
    cfl = 1;
    
    % iteration
    tt = 0:dt:t;
    xx = 0:dx:L;
    U = zeros(length(tt), length(xx));

    % first time step: the ghost node at i-1 is all zero since the bar was
    % at rest
    U(1,:) = 0; % initial condition
    U(2,:) = 0; % first step given that initial con is all zero
    U(2,1) = 1 - cos(2*dt);
    for ii = 3:length(tt)
        tnow = (ii-1)*dt;
        U(ii,:) = cfl*eval_uxx(U(ii-1,:),tnow,dt) + 2*U(ii-1,:) - U(ii-2,:);
        % boundary condition
        U(ii,end) = 0;
        if tnow <= pi
            U(ii,1) = 1 - cos(2*tnow);
        else
            U(ii,1) = 0;
        end
    end
    
end

function Uxx = eval_uxx(U, t, dt)
% evaluate d(du)/dxx
    % left boundary
    um1 = 1 - cos(2*(t-2*dt));
    up1 = 0;
    U = [um1, U, up1];
    Uxx = zeros(1,length(U)-2);
    for ii = 2:length(U)-1
        Uxx(ii-1) = U(ii-1) - 2*U(ii) + U(ii+1);
    end
end

