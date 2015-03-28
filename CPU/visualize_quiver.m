function [st] = visualize_quiver(index)
    [U, V] = read_velocity(index);
    [X, Y] = meshgrid(0:.01:.01*(size(U,1)-1), 0:.01:.01*(size(V,1)-1));
    st = quiver(X, Y, U, V);
    axis equal
    axis([X(1) X(end) Y(1) Y(end)]);
end
