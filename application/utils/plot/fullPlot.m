function fullPlot(y, y_est)
    figure
    plot3(y(1, :), y(2, :), y(3, :), 'b', ...
          y_est(1, :), y_est(2, :), y_est(3, :), 'rx', 'LineWidth', 2)
    xlabel('x (m)'), ylabel('y (m)'), zlabel('z (m)')
    title('Posição do Projétil')
    legend('Real', 'Medida', 'Filtrada', 'ItemHitFcn', @cb_legend)
    axis tight
    grid on
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1], 'Color', [1 1 1])
end