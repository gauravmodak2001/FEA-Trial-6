function elements_structure = create_elements_structure(data)
% PSEUDOCODE:
% 1. Extract data from input structure
% 2. Initialize variables and output structure
% 3. For each element:
%    a. Get nodal coordinates
%    b. Set up integration points
%    c. Calculate element stiffness matrix
%    d. Initialize surface force vector
%    e. Calculate shape functions and physical coordinates at integration points
%    f. Store element properties in output structure
% 4. Return element structures

connectivityMatrix = data.Connectivity;        % Element-node connectivity
allNodalLocations = data.GlobalNodes;          % Nodal coordinates
D = data.D;                                   % Material property matrix
nDOFPNode = data.nDOFPNode;                   % DOFs per node
n = 10;                                       % Grid size for visualization

numElements = size(connectivityMatrix, 1);     % Number of elements
numNodes = length(connectivityMatrix);         % Total number of nodes

% Initialize coordinate arrays
x(1:n+1, 1:n+1) = 0;                         % x-coordinates array
y(1:n+1, 1:n+1) = 0;                         % y-coordinates array

% Initialize output structure
elements_structure(numElements).shapeFunctions = [];
elements_structure(numElements).DiffShapefunctions = [];
elements_structure(numElements).Stiffnessmatrix = [];
elements_structure(numElements).SurfaceForce = [];
elements_structure(numElements).x = [];
elements_structure(numElements).y = [];

surface = 1;                                  % Default surface ID
% Commented out force node specifications
% Force_Nodes = [1];
% Force_Nodes = [0;0;0;0;0;0;0;1;1;1;1;1;1;1];

% Loop through each element
for i = 1:numElements
   % Get element data
   nNodes = size(connectivityMatrix, 2);      % Nodes per element
   nodeIndices = connectivityMatrix(i, :);    % Node indices for this element
   X = allNodalLocations(nodeIndices, :);     % Nodal coordinates
   X = reshape(X.', [], 1);                   % Reshape to column vector
   
   % Set up Gaussian quadrature points and weights
   Xi_coor = [0.93246, 0.661209, 0.23861, -0.23861, -0.661209, -0.93246];  % Integration points
   Weights = [0.1713244, 0.36076, 0.467913, 0.467913, 0.36076, 0.1713244]; % Integration weights
   
   % Alternative points (commented out)
   % Xi_coor = [-0.5077, 0.5077];
   % Weights = [1, 1];
   
   Gauss_in_Eta = 6;                          % Number of integration points in η direction
   Gauss_in_Xi = 6;                           % Number of integration points in ξ direction
   
   % Calculate element stiffness matrix: Kᵉ = ∫ BᵀDB dV
   Stiffnessmatrix = Calculate_K_matrix(nNodes, X, D, Gauss_in_Xi, Gauss_in_Eta, Xi_coor, Weights, data);
   
   % Initialize surface force vector
   % if(Force_Nodes(i,1) == 1)
   %     SurfaceForce = Surface_Traction_Calc(Gauss_in_Eta, nNodes, Xi_coor, X, Weights, surface, data);
   % else
   SurfaceForce = zeros(nDOFPNode * numNodes, 1);  % Zero surface force vector
   % end
   
   % Calculate shape functions and coordinates at integration points
   for j = 1:Gauss_in_Eta
       for k = 1:Gauss_in_Xi
           shapeFunctions = Shape_function([Xi_coor(j), Xi_coor(k)], nNodes);  % Shape functions
           DiffShapefunctions = Diff_Shape_functions([Xi_coor(j), Xi_coor(k)], nNodes);  % Shape function derivatives
           
           % Calculate physical coordinates at this point
           temp_C = calculate_coordinates(@Shape_function, [Xi_coor(j), Xi_coor(k)], X, nNodes);
           x(j, k) = temp_C(1, 1);  % x-coordinate
           y(j, k) = temp_C(2, 1);  % y-coordinate
       end
   end
   
   % Store element properties in output structure
   elements_structure(i).SurfaceForce = SurfaceForce;
   elements_structure(i).shapeFunctions = shapeFunctions;
   elements_structure(i).DiffShapefunctions = DiffShapefunctions;
   elements_structure(i).Stiffnessmatrix = Stiffnessmatrix;
   elements_structure(i).x = x;
   elements_structure(i).y = y;
end
end