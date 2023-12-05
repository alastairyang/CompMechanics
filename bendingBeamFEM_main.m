% solution from various polynomial power approximation
[x, w_5] = bendingBeamFEMFunc(5);
[~, w_2] = bendingBeamFEMFunc(2);
[~, w_3] = bendingBeamFEMFunc(3);
[~, w_4] = bendingBeamFEMFunc(4);
% define q: divide into two parts
q0 = 100; % N/m
L = 0.1; % m
x_1st = x(x<=L);
x_2nd = x(x<=2*L & x>L);
q_1st = q0*ones(size(x_1st));
q_fun = @(x)-(q0/L).*x + 2*q0; 
q_2nd = q_fun(x_2nd);
q = [q_1st, q_2nd];
% error in energy
err2 = calcEnergyNorm(w_5, w_2, x, q);
err3 = calcEnergyNorm(w_5, w_3, x, q);
err4 = calcEnergyNorm(w_5, w_4, x, q);
% plot
figure;
semilogy([2,3,4],[err2, err3, err4],'r-*','LineWidth',2)
xticks([2,3,4])
xlabel('Polynomial power','FontSize',15)
ylabel('Energy error','FontSize',15)
exportgraphics(gcf,'HW3_energError.png','Resolution',400)