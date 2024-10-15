function errorPlot(t_cell, errArray, plotTitle, plotLegend, xLim, yLim, lines)
    figure
    for i = 1:size(errArray, 1)
        subplot(1, size(errArray, 1), i)
        hold on
        for j = 1:size(errArray, 2)
            error = reshape(errArray{i, j}, 1, []);
            plot(t_cell{j, i}, error, lines(j), 'LineWidth', 2)
        end
        title(plotTitle(i), 'FontSize', 16)
        xlim(xLim)
        ylim(yLim)
        xlabel('tempo (s)', 'FontSize', 15), ylabel('erro (m)', 'FontSize', 12)
        lgd = legend(plotLegend, 'FontSize', 15, 'ItemHitFcn', @cb_legend);
        set(lgd, 'ItemTokenSize', [20, 50])
        grid on
    end
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1], 'Color', [1 1 1])
end