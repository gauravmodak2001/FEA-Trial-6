function u = Disp_function(Shape_function, Xi_coor, q, m)
% PSEUDOCODE:
% 1. Calculate shape functions at the given parametric coordinates
% 2. Use shape functions to interpolate displacements from nodal values
% 3. Return the displacement vector

% Function to calculate displacements using shape functions

% Calculate shape functions at given Xi coordinates
N = Shape_function(Xi_coor, m);            % Calculate shape function matrix N at point (ξ,η)

% Calculate displacements using the shape functions and input nodal displacements
% Formula: u = N·q where N contains shape functions and q contains nodal displacements
u = N * q;                                 % Interpolate displacements at point (ξ,η)
end