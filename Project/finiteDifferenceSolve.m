function [U, x] = finiteDifferenceSolve(xl, Ltot, N)
%FINITEDIFFERENCESOLVE

    h = Ltot/(N-1);
    x = xl + ([1:N]-1)*h;
    x = x(1:end-1);
    
    % Matrix [A]
    ee = ones(N,1);
    A = spdiags([ee -2*ee ee], -1:1, N, N);
    A = A(1:end-1,1:end-1);
    A(1,end) = 1; % periodic BC
    A(end,1) = 1; 

    % vector [b]
    func = @(x) -h^2*exp(-x.^2) + h^2*sqrt(pi)./(Ltot);
    b = transpose(func(x));
    
%     % solve
%     U = A\b;
    U = transpose(pinv(full(A))*b);

end

