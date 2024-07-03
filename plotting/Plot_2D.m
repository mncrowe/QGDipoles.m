function Plot_2D(F,x,y,plot_window,r)
% plots F(x,y)
%
% F: array with dimension 1 corresponding to x and dimension 2 to y
% (x,y): vectors
% plot_window: [min max] argument passed to xlim and ylim
% r: plot circle of radius r at origin

figure

pcolor(x,y,F')
shading interp
xlabel('x')
ylabel('y')

hold on

if nargin > 3
    xlim(plot_window)
    ylim(plot_window)

    line([plot_window(1) plot_window(1)],[plot_window(1) plot_window(end)],'Color','black')
    line([plot_window(end) plot_window(end)],[plot_window(1) plot_window(end)],'Color','black')
    line([plot_window(1) plot_window(end)],[plot_window(1) plot_window(1)],'Color','black')
    line([plot_window(1) plot_window(end)],[plot_window(end) plot_window(end)],'Color','black')
end

colormap(cmap(256,@(x) x.^1))
colorbar

contour(x,y,F',10,'k')

if nargin > 4
    fplot(@(t) r*sin(t), @(t) r*cos(t),[0,2*pi],'k--')
end

hold off

set(gca,'FontSize',12,'linewidth',0.7);

end