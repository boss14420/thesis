function [U, V] = read_velocity(index)
    U = dlmread(sprintf('uc%04d.mat', index));
    V = dlmread(sprintf('vc%04d.mat', index));
end
