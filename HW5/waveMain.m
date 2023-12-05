%% CEE6513 HW8
% Author: Donglai Yang
%% Sect 1: simulate
t = 100;
dt = 0.05;
L =1;

U = wave1DFunc(L,t,dt);

%% Sect 2: animate
dx = dt;
xx = 0:dx:L;

figure;
for i = 1:size(U,1)
    plot(xx, U(i,:),'k','LineWidth',2);
    ylim([-3,3])
    title(['t = ' num2str(i)])
    pause(0.01)
end

%% Sect 3: plot the 5 snapshots
t_plot = 20:20:100;
it = t_plot./dt;
dx = dt;
xx = 0:dx:L;
tt = 0:dt:t;

figure;
line_color = ["r","b","k","magenta","green"];
ct = 0;
for i = it
    ct = ct + 1;
    plot(xx, U(i,:),'Color',line_color(ct),'LineWidth',2); hold on;
    ylim([-3,3])
end
xlabel('Length','FontSize',14);
ylabel('Amplitude','FontSize',14)
legend(["20","40","60","80","100"])

exportgraphics(gcf,['HW5_plot_L' num2str(L) '.png'],'Resolution',300)