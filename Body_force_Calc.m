function [I] = Body_force_Calc(X, s, gauss_point_Xi, Integration_points, Weights, gauss_point_Eta, tz, f0, f1, f2, f3)
% PSEUDOCODE:
% 1. Initialize force vector
% 2. Loop through Gauss points in both directions:
%    a. Calculate shape functions at integration point
%    b. Calculate Jacobian and its determinant
%    c. Calculate physical coordinates at integration point
%    d. Evaluate body force vector
%    e. Add contribution to the force vector
% 3. Return body force vector

I = zeros(size(X));                        % Initialize body force vector

for i = 1:gauss_point_Xi                   % Loop over ξ direction points
   for j = 1:gauss_point_Eta              % Loop over η direction points
       % Calculate shape functions
       N = Shape_function([Integration_points(i), Integration_points(j)], s);
       
       % Calculate Jacobian and determinant
       J = Jacobian_Calculation(@Diff_Shape_functions, [Integration_points(i), Integration_points(j)], X, s);
       Det_J = det(J);                    % Determinant for area scaling
       
       % Calculate physical coordinates at this point
       X_vec = calculate_coordinates(@Shape_function, [Integration_points(i), Integration_points(j)], X, s);
       
       % Calculate body force vector (can vary with position)
       f_s = [f0 + (f1 * X_vec(1));       % x-component: f₀ + f₁·x
              f2 + (f3 * X_vec(2))];      % y-component: f₂ + f₃·y
       
       % Add contribution to force vector
       % Formula: I = ∫ N'·f |J| dξ dη
       I = I + N' * f_s * Det_J * Weights(i) * Weights(j) * tz;  % Weighted contribution
   end
end
end