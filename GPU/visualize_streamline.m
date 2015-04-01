function [st] = visualize_streamline(index)
    [U, V] = read_velocity(index);
    [X, Y] = meshgrid(0:.01:.01*(size(U,1)-1), 0:.01:.01*(size(V,1)-1));
    startX = zeros(1, size(U,1));
    startY = Y(:, 1);
    st = streamline(X, Y, U, V, startX, startY);
    axis equal
    axis([X(1) X(end) Y(1) Y(end)]);
end
