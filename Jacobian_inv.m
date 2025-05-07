function J_inv = Jacobian_inv(Diff_Shape_functions, Xi_coor, X, s)
% PSEUDOCODE:
% 1. Calculate Jacobian matrix J using Jacobian_Calculation
% 2. Determine number of dimensions from J
% 3. Calculate inverse of J
% 4. Create block diagonal matrix with J_Inverse repeated along diagonal
% 5. Return the block diagonal matrix

% Calculate Jacobian matrix J from parametric to physical coordinates:
% J = [∂x/∂ξ  ∂y/∂ξ]
%     [∂x/∂η  ∂y/∂η]
J = Jacobian_Calculation(Diff_Shape_functions, Xi_coor, X, s);

nDimensions = size(J, 1);             % Get number of dimensions (usually 2 for 2D problems)

% Calculate J⁻¹, the inverse of Jacobian matrix
% J⁻¹ = [∂ξ/∂x  ∂η/∂x]
%       [∂ξ/∂y  ∂η/∂y]
J_Inverse = inv(J);

% Create block diagonal matrix with repeated J_Inverse
% For 2D problems, resulting matrix is:
% J_inv = [∂ξ/∂x  ∂η/∂x  0      0    ]
%         [∂ξ/∂y  ∂η/∂y  0      0    ]
%         [0      0      ∂ξ/∂x  ∂η/∂x]
%         [0      0      ∂ξ/∂y  ∂η/∂y]
J_inv = repmat({J_Inverse}, 1, nDimensions);  % Create cell array with repeated J_Inverse
J_inv = blkdiag(J_inv{:});                    % Create block diagonal matrix
end