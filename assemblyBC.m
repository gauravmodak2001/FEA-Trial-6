function [K, F] = assemblyBC(data, Global_elements_structure)
% PSEUDOCODE:
% 1. Get maximum value from global stiffness matrix
% 2. Calculate penalty factor C
% 3. Extract global stiffness matrix and force vector
% 4. Loop through all boundary condition nodes:
%    a. Get DOF type and node number
%    b. Calculate global index
%    c. Apply penalty method to enforce boundary condition
% 5. Return modified stiffness matrix and force vector

% Find the maximum value in the stiffness matrix
k_max = max(Global_elements_structure.K, [], 'all');  % Maximum entry in K
C = k_max * 10^4;                                     % Penalty factor

% Extract global matrices
K = Global_elements_structure.K;                      % Global stiffness matrix
F = Global_elements_structure.F;                      % Global force vector

% Loop through all boundary conditions
for i = 1:length(data.BCnodes)
   nDOF = data.BCDof(i);                            % DOF type (1=x, 2=y)
   nNode = data.BCnodes(i);                         % Node number
   
   % Calculate global index in system
   index = convert_to_index(nNode, nDOF, data.nDOFPNode);
   
   % Apply boundary conditions using penalty method
   % Modifies equation K·u = F to enforce u = value at constrained DOFs
   F(index) = F(index) + C * data.BC(i);            % Add large force = C·value
   K(index, index) = K(index, index) + C;           % Add large stiffness
end
end