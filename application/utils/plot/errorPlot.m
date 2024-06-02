function errorPlot(t, errArray, plotTitle, plotLegend, ylim_max)
    figure
    for i = 1:size(errArray, 1)
        subplot(1, size(errArray, 1), i)
        hold on
        for j = 1:size(errArray, 2)
            error = reshape(errArray(i, j, :), 1, []);
            plot(t, error, 'LineWidth', 2)
        end
        title(plotTitle(i))
        ylim([0 ylim_max])
        xlabel('tempo (s)'), ylabel('erro (m)')
        legend(plotLegend, 'ItemHitFcn', @cb_legend)
        grid on
    end
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1], 'Color', [1 1 1]);
end