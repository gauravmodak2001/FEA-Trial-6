function dN = Diff_Shape_functions(Xi_Coor, s)
% PSEUDOCODE:
% 1. Check element type (4-node or 8-node)
% 2. If 4-node element:
%    a. Calculate derivatives of 4 shape functions
%    b. Arrange derivatives in matrix format
% 3. If 8-node element:
%    a. Calculate derivatives of 8 shape functions
%    b. Arrange derivatives in matrix format
% 4. Return matrix of shape function derivatives

if s == 4
   % 4-node configuration for the position (-1,1)
   xi = Xi_Coor(1, 1);                      % Extract ξ coordinate
   eta = Xi_Coor(1, 2);                     % Extract η coordinate
   
   % Calculate derivatives ∂Nᵢ/∂ξ and ∂Nᵢ/∂η for (N1,N2,N3,N4)
   dN1dxi = eta / 4 - 1 / 4;                % ∂N₁/∂ξ = (η-1)/4
   dN1deta = xi / 4 - 1 / 4;                % ∂N₁/∂η = (ξ-1)/4
   
   dN2dxi = 1 / 4 - eta / 4;                % ∂N₂/∂ξ = (1-η)/4
   dN2deta = -xi / 4 - 1 / 4;               % ∂N₂/∂η = (-ξ-1)/4
   
   dN3dxi = eta / 4 + 1 / 4;                % ∂N₃/∂ξ = (η+1)/4
   dN3deta = xi / 4 + 1 / 4;                % ∂N₃/∂η = (ξ+1)/4
   
   dN4dxi = -eta / 4 - 1 / 4;               % ∂N₄/∂ξ = (-η-1)/4
   dN4deta = 1 / 4 - xi / 4;                % ∂N₄/∂η = (1-ξ)/4
   
   % Combine derivatives into dN matrix for strain calculations
   dN = [dN1dxi, 0, dN2dxi, 0, dN3dxi, 0, dN4dxi, 0;
         dN1deta, 0, dN2deta, 0, dN3deta, 0, dN4deta, 0;
         0, dN1dxi, 0, dN2dxi, 0, dN3dxi, 0, dN4dxi;
         0, dN1deta, 0, dN2deta, 0, dN3deta, 0, dN4deta];
         
elseif s == 8
   % 8-node configuration for the position (-1,1)
   Xi = Xi_Coor(1, 1);                      % Extract ξ coordinate
   Eta = Xi_Coor(1, 2);                     % Extract η coordinate
   
   % Calculate derivatives for 8-node serendipity element
   dN1dxi = -((Eta + 2 * Xi) * (Eta - 1)) / 4;    % ∂N₁/∂ξ
   dN1deta = -((2 * Eta + Xi) * (Xi - 1)) / 4;    % ∂N₁/∂η
   
   dN2dxi = Xi * (Eta - 1);                       % ∂N₂/∂ξ
   dN2deta = ((Xi - 1) * (Xi + 1)) / 2;           % ∂N₂/∂η
   
   dN3dxi = ((Eta - 2 * Xi) * (Eta - 1)) / 4;     % ∂N₃/∂ξ
   dN3deta = ((Xi + 1) * (2 * Eta - Xi)) / 4;     % ∂N₃/∂η
   
   dN4dxi = -((Eta - 1) * (Eta + 1)) / 2;         % ∂N₄/∂ξ
   dN4deta = -Eta * (Xi + 1);                     % ∂N₄/∂η
   
   dN5dxi = ((Eta + 2 * Xi) * (Eta + 1)) / 4;     % ∂N₅/∂ξ
   dN5deta = ((2 * Eta + Xi) * (Xi + 1)) / 4;     % ∂N₅/∂η
   
   dN6dxi = -Xi * (Eta + 1);                      % ∂N₆/∂ξ
   dN6deta = -((Xi - 1) * (Xi + 1)) / 2;          % ∂N₆/∂η
   
   dN7dxi = -((Eta - 2 * Xi) * (Eta + 1)) / 4;    % ∂N₇/∂ξ
   dN7deta = -((Xi - 1) * (2 * Eta - Xi)) / 4;    % ∂N₇/∂η
   
   dN8dxi = ((Eta - 1) * (Eta + 1)) / 2;          % ∂N₈/∂ξ
   dN8deta = Eta * (Xi - 1);                      % ∂N₈/∂η
   
   % Combine derivatives into dN matrix for strain calculations
   dN = [dN1dxi, 0, dN2dxi, 0, dN3dxi, 0, dN4dxi, 0, dN5dxi, 0, dN6dxi, 0, dN7dxi, 0, dN8dxi, 0;
         dN1deta, 0, dN2deta, 0, dN3deta, 0, dN4deta, 0, dN5deta, 0, dN6deta, 0, dN7deta, 0, dN8deta, 0;
         0, dN1dxi, 0, dN2dxi, 0, dN3dxi, 0, dN4dxi, 0, dN5dxi, 0, dN6dxi, 0, dN7dxi, 0, dN8dxi;
         0, dN1deta, 0, dN2deta, 0, dN3deta, 0, dN4deta, 0, dN5deta, 0, dN6deta, 0, dN7deta, 0, dN8deta];
else
   error('Invalid input. Please enter either 4 or 8 for the number of nodes.');
end
end