function Displacementandstraincalc = Displacementandstraincalc(Global_elements_structure, data)
% PSEUDOCODE:
% 1. Setup mesh in parametric space for result visualization
% 2. Initialize arrays for strain, stress, and displacements
% 3. For each element:
%    a. Get element connectivity and nodal coordinates
%    b. Extract element displacement vector from global solution
%    c. For each visualization point in parametric space:
%       i.   Calculate B matrix at this point
%       ii.  Calculate strain: ε = B·q
%       iii. Calculate stress: σ = D·ε
%       iv.  Calculate physical coordinates (x,y)
%       v.   Calculate displacements: u = N·q
%       vi.  Store results in output structure
% 4. Return calculated results

connectivityMatrix = data.Connectivity;           % Get element-node connectivity
allNodalLocations = data.updateNodalLocations;    % Get nodal coordinates

% Define parametric space grid for visualization
x_min = -1;                                       % Lower bound of parametric space
x_max = 1;                                        % Upper bound of parametric space
n = 10;                                           % Number of divisions in each direction

% Initialize arrays for parametric coordinates
Xi = zeros(n + 1, 1);                            % Xi coordinates array
Eta = zeros(n + 1, 1);                           % Eta coordinates array

% Initialize arrays for results
Strain = zeros(n + 1, n + 1, 3);                 % Strain tensor array (εxx, εyy, γxy)
Stress = zeros(n + 1, n + 1, 3);                 % Stress tensor array (σxx, σyy, τxy)
x(1:n+1, 1:n+1) = 0;                             % Physical x-coordinates
y(1:n+1, 1:n+1) = 0;                             % Physical y-coordinates
Displacements = zeros(n + 1, n + 1, 2);          % Displacement vector array (u, v)

% Get element and node counts
numElements = size(connectivityMatrix, 1);        % Number of elements
numNodes = size(connectivityMatrix, 2);           % Number of nodes per element

% Initialize output structure
Displacementandstraincalc(numElements).U = [];
Displacementandstraincalc(numElements).Strain = [];
Displacementandstraincalc(numElements).Stress = [];
Displacementandstraincalc(numElements).Xi = [];
Displacementandstraincalc(numElements).Eta = [];
Displacementandstraincalc(numElements).x = [];
Displacementandstraincalc(numElements).y = [];

% Loop through each element
for i = 1:numElements
   % Get element nodal data
   nodeIndices = connectivityMatrix(i, :);                  % Global node indices for this element
   X = allNodalLocations(nodeIndices, :);                   % Nodal coordinates
   X = reshape(X.', [], 1);                                 % Reshape to column vector [x1;y1;x2;y2;...]
   
   % Extract element displacement vector from global solution
   % Formula: q_e = T_e·q where T_e is the transformation matrix
   Q_indices = (nodeIndices - 1) * 2 + 1;                   % Starting DOF indices for each node
   Q = Global_elements_structure.Q([Q_indices; Q_indices + 1]); % Extract displacements [u1;v1;u2;v2;...]
   Q = reshape(Q, [], 1);                                   % Reshape to column vector
   
   % Loop through visualization grid in parametric space
   for j = 1:(n + 1)
       Xi(j, 1) = x_min(1) + ((j - 1) * (x_max(1) - x_min(1))) / n;  % Calculate ξ coordinate
       
       for m = 1:(n + 1)
           Eta(m, 1) = x_min(1) + ((m - 1) * (x_max(1) - x_min(1))) / n;  % Calculate η coordinate
           
           % Calculate B matrix at this parametric point
           B = Calculate_B_matrix(@Diff_Shape_functions, [Xi(j, 1), Eta(m, 1)], X, numNodes);
           
           % Calculate strain: ε = B·q
           temp_V = strain_cal(B, Q);
           
           % Store strain and stress: σ = D·ε
           Strain(j, m, :) = temp_V;
           Stress(j, m, :) = data.D * temp_V;
           
           % Calculate physical coordinates: x = N·X
           temp_C = calculate_coordinates(@Shape_function, [Xi(j, 1), Eta(m, 1)], X, numNodes);
           x(j, m) = temp_C(1, 1);
           y(j, m) = temp_C(2, 1);
           
           % Calculate displacements: u = N·q
           temp_U = Disp_function(@Shape_function, [Xi(j, 1), Eta(m, 1)], Q, numNodes);
           
           % Store displacement components
           for o = 1:2
               Displacements(j, m, o) = temp_U(o, 1);
               Displacementandstraincalc(i).U(j, m, o) = Displacements(j, m, o);
           end
           
           % Store strain and stress vectors
           strain_vector = reshape(Strain(j, m, :), [], 1);
           Stress_vector = reshape(Stress(j, m, :), [], 1);
           Displacementandstraincalc(i).Strain(j, m, :) = strain_vector;
           Displacementandstraincalc(i).Stress(j, m, :) = Stress_vector;
       end
   end
   
   % Store parametric and physical coordinates for this element
   Displacementandstraincalc(i).Xi = Xi;
   Displacementandstraincalc(i).Eta = Eta;
   Displacementandstraincalc(i).x = x;
   Displacementandstraincalc(i).y = y;
end
end