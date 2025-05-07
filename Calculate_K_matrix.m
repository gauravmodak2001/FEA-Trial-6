function [K] = Calculate_K_matrix(nNodes, X, D, gauss_point_Xi, gauss_point_Eta, Integration_points, Weights, data)
% PSEUDOCODE:
% 1. Initialize matrices for integration and stiffness
% 2. Loop through Gauss points in both directions:
%    a. Calculate B matrix at each integration point
%    b. Calculate Jacobian and its determinant
%    c. Evaluate integrand (B'DB) and multiply by weight and determinant
%    d. Accumulate contributions to stiffness matrix
% 3. Multiply by thickness to get final stiffness matrix
% 4. Return element stiffness matrix

% Initialize I matrix for integration
I = zeros(size(X, 1));                       % Initialize integration accumulator
K = zeros(size(X, 1));                       % Initialize stiffness matrix

% Double loop over integration points
for i = 1:gauss_point_Xi                     % Loop over ξ direction points
   for j = 1:gauss_point_Eta                % Loop over η direction points
       % Calculate B matrix: relates displacement to strain, ε = B·u
       B = Calculate_B_matrix(@Diff_Shape_functions, [Integration_points(i), Integration_points(j)], X, nNodes);
       
       % Calculate Jacobian: maps from parametric to physical space
       J = Jacobian_Calculation(@Diff_Shape_functions, [Integration_points(i), Integration_points(j)], X, nNodes);
       
       % Calculate determinant of J: scaling factor for area/volume
       Det_J = det(J);
       
       % Evaluate function: B'DB represents material stiffness
       % Formula: Kᵉ = ∫ B'DB |J| dξ dη
       F = Weights(i) * Weights(j) * B' * D * B * Det_J;  % Weighted contribution
       I = I + F;                           % Accumulate contributions
       K = I * data.thikness;               % Multiply by thickness for plane problems
   end
end
end