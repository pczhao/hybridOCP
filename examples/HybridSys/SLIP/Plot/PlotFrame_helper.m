
if ~exist('mycolor','var')
    mycolor = [ 0, 0.4470, 0.7410 ];
end
h_2d = gca;

hold on;
idx_mat = mode_hist(2:end) - mode_hist(1:end-1);
idx = [ find( idx_mat > 0 ); find( idx_mat < 0 )+1];
idx = [1; idx; length(mode_hist)-1];
springcoord( [ 0, 0 ] , [ 0, params.l0 ], 5, 0.3, 0.04);

for i = 1 : length(idx)
    j = idx(i);
    if mode_hist(j+1) == mode_hist(j)
        xv = (x_hist(j+1) - x_hist(j)) / (t_hist(j+1) - t_hist(j));
        yv = (y_hist(j+1) - y_hist(j)) / (t_hist(j+1) - t_hist(j));
    elseif mode_hist(j) == mode_hist(j-1)
        xv = (x_hist(j) - x_hist(j-1)) / (t_hist(j) - t_hist(j-1));
        yv = (y_hist(j) - y_hist(j-1)) / (t_hist(j) - t_hist(j-1));
    end
    PlotFrame(state_hist(j,:),params,xv,yv,h_2d);
end

plot( h_2d, x_hist, y_hist, '-', 'color', mycolor, 'LineWidth', 1);

axis equal
ylim([0, 0.35]);
xlim([-1,0]);
xlabel('$a$','Interpreter','LaTex','FontSize',15);
ylabel('$b$','Interpreter','LaTex','FontSize',15);
box on;
set(h_2d, 'YTick', [0 1.5]);