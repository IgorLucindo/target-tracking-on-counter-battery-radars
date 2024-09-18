function errorPlotStd(t_cell, errArray, errArrayStd, plotTitle, xLim, yLim_max, lines)
    figure
    for i = 1:size(errArray, 1)
        subplot(1, size(errArray, 1), i)
        hold on
        for j = 1:size(errArray, 2)
            error = reshape(errArray{i, j}, 1, []);
            fill([t_cell{j, i}, fliplr(t_cell{j, i})], ...
                [error + errArrayStd{i, j}, ...
                fliplr(error - errArrayStd{i, j})], ...
                [0.8 0.8 1], 'EdgeColor', [0.7 0.7 1]);
            hold on
            plot(t_cell{j, i}, error, lines(j), 'LineWidth', 2)
        end
        title(plotTitle(i))
        xlim(xLim)
        ylim([0 yLim_max])
        xlabel('tempo (s)'), ylabel('erro (m)')
        grid on
    end
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1], 'Color', [1 1 1])
end