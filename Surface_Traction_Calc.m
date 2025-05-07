function [I] = Surface_Traction_Calc(gauss_point_Xi, s, Integration_points, X, Weights, surface, data)
% Surface_Traction_Calc - Calculates surface traction forces for FEM analysis
%
% PSEUDOCODE:
% 1. Initialize force vector I with zeros
% 2. Loop through Gauss integration points
% 3. For each integration point:
%    a. Based on the surface number (1-4), determine parametric coordinates
%    b. Calculate shape functions N at those coordinates
%    c. Calculate Jacobian J at those coordinates
%    d. Calculate determinant of J for the surface segment
%    e. Calculate physical coordinates X_vec at those parametric coordinates
%    f. Compute force vector f_s based on force parameters and coordinates
%    g. Add contribution to force vector: N'*f_s*Det_J*Weight*thickness
% 4. Return the final surface force vector I
%
% Input:
%   gauss_point_Xi - Number of Gauss integration points
%   s - Shape function parameter (number of nodes in element)
%   Integration_points - Array of Gauss quadrature points
%   X - Nodal coordinates
%   Weights - Integration weights for Gauss quadrature
%   surface - Surface number (1=right, 2=left, 3=top, 4=bottom)
%   data - Structure containing force and thickness information
%
% Output:
%   I - Calculated surface traction force vector

I = zeros;                            % Initialize force vector I with zeros

for i = 1:gauss_point_Xi              % Loop through Gauss integration points
    if surface == 1                   % Surface 1 (right edge: Xi=1)
        N = Shape_function([1, Integration_points(i)], s);  % Calculate shape functions at parametric coordinates
        J = Jacobian_Calculation(@Diff_Shape_functions, [1, Integration_points(i)], X, s);  % Calculate Jacobian
        Det_J = sqrt(J(2,1)^2 + J(2,2)^2);  % Calculate determinant for length of surface segment
        X_vec = calculate_coordinates(@Shape_function, [1, Integration_points(i)], X, s);  % Get physical coordinates
        f_s = [data.f(1) + data.f(2) * X_vec(1); data.f(3)*X_vec(2) + data.f(4)*X_vec(2)];  % Calculate force vector
        % Alternative calculation (commented out):
        % f_s = [data.f(1) + (data.f(2) * X_vec(1)); data.f(3) + (data.f(4)*X_vec(2))];
        I = I + N' * f_s * Det_J * Weights(i) * data.thikness;  % Add contribution to force vector
        
    elseif surface == 2               % Surface 2 (left edge: Xi=-1)
        N = Shape_function([-1, Integration_points(i)], s);  % Calculate shape functions
        J = Jacobian_Calculation(@Diff_Shape_functions, [-1, Integration_points(i)], X, s);  % Calculate Jacobian
        Det_J = sqrt(J(2,1)^2 + J(2,2)^2);  % Calculate determinant
        X_vec = calculate_coordinates(@Shape_function, [-1, Integration_points(i)], X, s);  % Get physical coordinates
        f_s = [data.f(1) + (data.f(2) * X_vec(1)); data.f(3) + (data.f(4)*X_vec(2))];  % Calculate force vector
        I = I + N' * f_s * Det_J * Weights(i) * data.thikness;  % Add contribution
        
    elseif surface == 3               % Surface 3 (top edge: Eta=1)
        % Note: commented indicator with empty comment
        N = Shape_function([Integration_points(i), 1], s);  % Calculate shape functions
        J = Jacobian_Calculation(@Diff_Shape_functions, [Integration_points(i), 1], X, s);  % Calculate Jacobian
        Det_J = sqrt(J(1,1)^2 + J(1,2)^2);  % Calculate determinant
        X_vec = calculate_coordinates(@Shape_function, [Integration_points(i), 1], X, s);  % Get physical coordinates
        f_s = [data.f(1) + data.f(2) * X_vec(1); data.f(3)*X_vec(2) + data.f(4)*X_vec(1)];  % Calculate force vector
        % Alternative calculation (commented out):
        % f_s = [data.f(1) + (data.f(2) * X_vec(1)); data.f(3) + (data.f(4)*X_vec(2))];
        I = I + N' * f_s * Det_J * Weights(i) * data.thikness;  % Add contribution
        
    elseif surface == 4               % Surface 4 (bottom edge: Eta=-1)
        N = Shape_function([Integration_points(i), -1], s);  % Calculate shape functions
        J = Jacobian_Calculation(@Diff_Shape_functions, [Integration_points(i), -1], X, s);  % Calculate Jacobian
        Det_J = sqrt(J(1,1)^2 + J(1,2)^2);  % Calculate determinant
        X_vec = calculate_coordinates(@Shape_function, [Integration_points(i), -1], X, s);  % Get physical coordinates
        f_s = [data.f(1) + (data.f(2) * X_vec(1)); data.f(3) + (data.f(4)*X_vec(2))];  % Calculate force vector
        I = I + N' * f_s * Det_J * Weights(i) * data.thikness;  % Add contribution
    end
end
end