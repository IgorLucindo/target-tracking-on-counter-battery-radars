function dynamicPlotRecord(y_true, y_est, impactPointPredArray, shootingPointPredArray)
    video = VideoWriter('plot.mp4', 'MPEG-4');
    video.FrameRate = 60;
    open(video)

    zero_columns = all(impactPointPredArray == 0, 1);
    num_zero_columns = sum(zero_columns);

    figure

    plot(y_true(2, :), y_true(3, :), 'b', 'LineWidth', 2)
    hold on

    h = animatedline('LineStyle', '-', 'Color', 'r', 'LineWidth', 4);
    p1 = animatedline('Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
    p2 = animatedline('Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 5);

    axis([-1e3, 12e3, 0, 4e3])
    %view(45, 15);
    xlabel('x (m)'), ylabel('z (m)'), zlabel('z (m)')
    title('Posição do Projétil')
    legend('Medida', 'Filtrada', 'Ponto de Impacto', 'Ponto de Disparo')
    grid on
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1], 'Color', [1 1 1]);

    for i = 1:num_zero_columns
        addpoints(h, y_est(2, i), y_est(3, i));
    end
    for i = num_zero_columns + 1:length(impactPointPredArray)
        clearpoints(p1);
        clearpoints(p2);
        addpoints(h, y_est(2, i), y_est(3, i));
        addpoints(p1, impactPointPredArray(2, i), impactPointPredArray(3, i));
        addpoints(p2, shootingPointPredArray(2, i), shootingPointPredArray(3, i));
        drawnow
        frame = getframe(gcf);
        writeVideo(video, frame)
    end

    close(video)
end