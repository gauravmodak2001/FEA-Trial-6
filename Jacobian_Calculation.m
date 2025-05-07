function J = Jacobian_Calculation(Diff_Shape_functions, Xi_Coor, X, s)
% PSEUDOCODE:
% 1. Calculate derivatives of shape functions at parametric coordinates
% 2. Compute Jacobian matrix by multiplying shape function derivatives with nodal coordinates
% 3. Reshape result to proper square matrix form
% 4. Return the Jacobian matrix

% Calculate the derivative of shape functions dN/dξ and dN/dη
% dN_cap contains [∂N₁/∂ξ, ∂N₂/∂ξ, ...; ∂N₁/∂η, ∂N₂/∂η, ...]
dN_cap = Diff_Shape_functions(Xi_Coor, s);

% Compute the Jacobian matrix: J = [∂x/∂ξ ∂y/∂ξ; ∂x/∂η ∂y/∂η]
% J = dN_cap * X where X contains nodal coordinates
J = dN_cap * X;

% Reshape the Jacobian matrix to a square matrix
matrixSize = sqrt(size(J, 1));
J = reshape(J, [matrixSize, matrixSize]);
end