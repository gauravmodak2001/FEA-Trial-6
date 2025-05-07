% FEM Analysis for 8-noded 2D element assembly
% Problem parameters: E=70,000, v=0.33, plane stress, t_z = 0.3
% Surface force fs = [0; 15] at the top edge

% Clear workspace
clear all; close all; clc;

% Set up data structure
data = struct();

% Define material properties
data.E = 70000;        % Young's modulus
data.v = 0.33;         % Poisson's ratio  
data.thikness = 0.3;   % Thickness
data.nDOFPNode = 2;    % Degrees of freedom per node (x,y)

% Define plane stress condition
data.D = D_matrix_Stress(data);

% Define node coordinates (26 nodes total)
data.GlobalNodes = [
    0, 6;    % Node 1
    1.5, 6;  % Node 2
    3, 6;    % Node 3
    4.5, 6;  % Node 4
    6, 6;    % Node 5
    0, 4.5;  % Node 6
    3, 4.5;  % Node 7
    6, 4.5;  % Node 8
    0, 3;    % Node 9
    1.5, 3;  % Node 10
    3, 3;    % Node 11
    4.5, 3;  % Node 12
    6, 3;    % Node 13
    0, 1.5;  % Node 14
    3, 1.5;  % Node 15
    0, 0;    % Node 16
    1.5, 0;  % Node 17
    3, 0;    % Node 18
    3.5, 0;  % Node 19
    3.792893, 2.207107; % Node 20
    6, 2.5;  % Node 21
    4, 0;    % Node 22
    4.152241, 0.765367; % Node 23
    4.585786, 1.414214; % Node 24
    5.234633, 1.847759; % Node 25
    6, 2;    % Node 26
];

% Element connectivity with 8-node elements (5 elements total)
data.Connectivity = [
    9, 10, 11, 7, 3, 2, 1, 6;      % Element 1
    11, 12, 13, 8, 5, 4, 3, 7;      % Element 2
    16, 17, 18, 15, 11, 10, 9, 14;  % Element 3
    18, 19, 22, 23, 24, 20, 11, 15; % Element 4
    24, 25, 26, 21, 13, 12, 11, 20; % Element 5
];

% Set the nodes per element
data.numElements = size(data.Connectivity, 1);
data.nTotalDOF = size(data.GlobalNodes, 1) * data.nDOFPNode;
data.type_of_mesh = "Undistorted";
data.updateNodalLocations = data.GlobalNodes;

% Create a figure to visualize the mesh
figure;
hold on;
title('Element Connectivity Visualization');
xlabel('X-coordinate');
ylabel('Y-coordinate');
grid on;

% Define colors for each element
colors = ['r', 'g', 'b', 'c', 'm'];

% Plot each element with a different color
for i = 1:data.numElements
    elementNodes = data.Connectivity(i, :);
    nodeCoords = data.GlobalNodes(elementNodes, :);
    
    % Close the element by repeating the first node
    nodeCoords = [nodeCoords; nodeCoords(1, :)];
    
    % Plot element with specific color
    plot(nodeCoords(:, 1), nodeCoords(:, 2), colors(mod(i-1, length(colors))+1), 'LineWidth', 2);
    
    % Label the nodes
    for j = 1:length(elementNodes)
        text(data.GlobalNodes(elementNodes(j), 1), data.GlobalNodes(elementNodes(j), 2), ...
            num2str(elementNodes(j)), 'FontSize', 10, 'Color', 'b');
    end
end

% Set boundary conditions
% Nodes 16, 17, 18, 19, 22 are fixed (bottom nodes)
% Nodes 5, 8, 13, 21, 26 are fixed (right side nodes)
fixedNodesBottom = [16, 17, 18, 19, 22];
fixedNodesRight = [5, 8, 13, 21, 26];
fixedNodes = [fixedNodesBottom, fixedNodesRight];
numFixedNodes = length(fixedNodes);

% Initialize BC arrays
data.BCnodes = zeros(numFixedNodes * 2, 1);
data.BCDof = zeros(numFixedNodes * 2, 1);
data.BC = zeros(numFixedNodes * 2, 1);

% Set up fixed boundary conditions (both x and y directions)
for i = 1:numFixedNodes
    % X-direction constraint
    data.BCnodes(2*i-1) = fixedNodes(i);
    data.BCDof(2*i-1) = 1;
    data.BC(2*i-1) = 0;
    
    % Y-direction constraint
    data.BCnodes(2*i) = fixedNodes(i);
    data.BCDof(2*i) = 2;
    data.BC(2*i) = 0;
end

% Set up surface force parameters
data.f = [0, 0, 0, 0]; % No body forces

% Create elements structure
elements_structure = create_elements_structure(data);

% Set up surface force parameters for top edge (nodes 1-5)
data.f = [0, 0, 15, 0]; % Surface force [0; 15] in y-direction

% Create an array to track which elements have surface loads
elementsWithSurfaceLoad = zeros(data.numElements, 1);

% Elements 1 and 2 have the top edge (where force is applied)
elementsWithSurfaceLoad(1) = 1; % Element 1 has top edge with nodes 1,2,3
elementsWithSurfaceLoad(2) = 1; % Element 2 has top edge with nodes 3,4,5

% Now update the element structure with surface forces
for i = 1:data.numElements
    if elementsWithSurfaceLoad(i) == 1
        % Get the element nodes
        elementNodes = data.Connectivity(i, :);
        % Get coordinates for this element
        X = reshape(data.GlobalNodes(elementNodes, :)', [], 1);
        
        % Calculate surface force for the top edge (surface 3 corresponds to top edge in our convention)
        % Use Integration_points and Weights for Gauss quadrature
        [Integration_points, Weights] = Gauss_quadrant(3); % 3 Gauss points for accuracy
        
        % Calculate surface traction
        surface = 3; % Top edge
        surfaceForce = Surface_Traction_Calc(3, 8, Integration_points, X, Weights, surface, data);
        
        % Add this to the element's surface force
        elements_structure(i).SurfaceForce = elements_structure(i).SurfaceForce + surfaceForce;
    end
end

% Create global stiffness matrix and force vector
Global_elements_structure = Global_Stifness_Matrix(@convert_to_index, elements_structure, data);

% Apply boundary conditions
[K, F] = assemblyBC(data, Global_elements_structure);

% Solve the system (K*Q = F)
Q = K \ F;
Global_elements_structure.Q = Q;

% Display the global stiffness matrix and force vector
disp('Global Stiffness Matrix K:');
disp(K);

disp('Global Force Vector F:');
disp(F);

% Create a heatmap of the stiffness matrix
figure('Position', [100, 100, 800, 700]);
imagesc(log10(abs(K) + 1e-10)); % Use log scale for better visualization
colormap(jet); % Use jet colormap for better contrast
colorbar; % Add colorbar
title('Heat Map of Global Stiffness Matrix K', 'FontSize', 14);
xlabel('DOF Index', 'FontSize', 12);
ylabel('DOF Index', 'FontSize', 12);
axis equal tight;
set(gca, 'FontSize', 10);

% Add grid lines
hold on;
[rows, cols] = size(K);
for i = 0.5:1:rows+0.5
    line([0.5, cols+0.5], [i, i], 'Color', [0.5, 0.5, 0.5, 0.2], 'LineWidth', 0.5);
end
for j = 0.5:1:cols+0.5
    line([j, j], [0.5, rows+0.5], 'Color', [0.5, 0.5, 0.5, 0.2], 'LineWidth', 0.5);
end

% Add text for the colorbar
cb = colorbar;

% Find non-zero components of the force vector
nonZeroIndices = find(abs(F) > 1e-10);
fprintf('Non-zero components of the global force vector F:\n');
for i = 1:length(nonZeroIndices)
    fprintf('F(%d) = %.6f\n', nonZeroIndices(i), F(nonZeroIndices(i)));
end

% Plot the force vector
figure('Position', [100, 100, 800, 400]);
stem(F, 'filled', 'LineWidth', 1.5);
grid on;
title('Global Force Vector F', 'FontSize', 14);
xlabel('DOF Index', 'FontSize', 12);
ylabel('Force Value', 'FontSize', 12);
set(gca, 'FontSize', 10);

% Add text labels for significant force values
threshold = max(abs(F))/100; % Threshold to show labels
for i = 1:length(F)
    if abs(F(i)) > threshold
        text(i, F(i), sprintf('%.1f', F(i)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', 8);
    end
end

% Highlight non-zero force components
hold on;
stem(nonZeroIndices, F(nonZeroIndices), 'r', 'filled', 'LineWidth', 2);
hold off;

% Add legend
legend('All DOF Forces', 'Non-zero Forces', 'Location', 'best');

