function [B_matrix_result] = Calculate_B_matrix(Diff_Shape_functions, Xi_Coor, X, s)
% PSEUDOCODE:
% 1. Define A matrix to extract physical strain components
% 2. Calculate derivatives of shape functions in parametric space
% 3. Calculate inverse Jacobian to transform derivatives
% 4. Compute B matrix by combining these components
% 5. Return B matrix for strain calculations

% Define matrix A_matrix for strain-displacement relationship
A_matrix = [1 0 0 0;                    % Component for εxx
            0 0 0 1;                     % Component for εyy
            0 1 1 0];                    % Component for γxy (engineering shear strain)

% Compute the derivative of shape functions
Derivative_Nhat_Dxi = Diff_Shape_functions(Xi_Coor, s);  % ∂N/∂ξ and ∂N/∂η

% Compute the inverse Jacobian matrix to transform derivatives
J_hat_Inverse = Jacobian_inv(@Diff_Shape_functions, Xi_Coor, X, s);  % J⁻¹

% Compute B matrix: B = A · J⁻¹ · ∂N/∂ξ
% This maps from displacements to strains: ε = B·q
B_matrix_result = A_matrix * J_hat_Inverse * Derivative_Nhat_Dxi;
end