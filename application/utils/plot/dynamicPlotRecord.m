function dynamicPlotRecord(y1, y2)
    video = VideoWriter('plot.mp4', 'MPEG-4');
    video.FrameRate = 60;
    open(video)

    figure

    h1 = animatedline('LineStyle', '-', 'Color', 'r', 'LineWidth', 4);
    h2 = animatedline('LineStyle', '-', 'Color', 'b', 'LineWidth', 4);

    axis([-1e3, 12e3, -1e3, 12e3, 0, 4e3])
    view(45, 15);
    xlabel('x (m)'), ylabel('z (m)'), zlabel('z (m)')
    title('Posição do Projétil')
    legend("trajetória parabólica (modelo linear)", "trajetória com arrasto (modelo não linear)")
    grid on
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1], 'Color', [1 1 1]);

    for i = 1:length(y1)
        addpoints(h1, y1(1, i), y1(2, i), y1(3, i));
        drawnow
        if i <= length(y2)
            addpoints(h2, y2(1, i), y2(2, i), y2(3, i));
            drawnow
        end
        frame = getframe(gcf);
        writeVideo(video, frame)
    end

    close(video)
end