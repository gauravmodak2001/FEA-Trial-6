function [x] = calculate_coordinates(Shape_function, Xi_Coor, X, s)
% PSEUDOCODE:
% 1. Calculate shape functions at given parametric coordinates
% 2. Interpolate physical coordinates using shape functions
% 3. Return physical coordinates

% Calling the shape function
N = Shape_function(Xi_Coor, s);           % Calculate shape functions at parametric point (ξ,η)

% Calculating the physical coordinates
% Formula: x = N·X maps from parametric to physical space
x = N * X;                               % Interpolate physical coordinates
end