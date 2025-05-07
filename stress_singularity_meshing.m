% Stress Singularity Analysis with Multiple Mesh Refinements
% Clear workspace
clear all; close all; clc;

% Define material properties
E = 70000;        % Young's modulus
v = 0.33;         % Poisson's ratio  
thickness = 0.3;  % Thickness
problemWidth = 2; % Width of rectangle (x-direction)
problemHeight = 6; % Height of rectangle (y-direction)

% Array of mesh sizes (each one half the previous)
meshSizes = [1, 0.5, 0.25, 0.125];

% Create a cell array to store results
modelResults = cell(length(meshSizes), 1);

% Loop through each mesh size and create models
for i = 1:length(meshSizes)
    % Create model with current mesh size
    fprintf('\n===== Processing mesh size = %f =====\n', meshSizes(i));
    
    % Initialize data structure (similar to original code)
    data = struct();
    data.E = E;
    data.v = v;
    data.thikness = thickness;
    data.nDOFPNode = 2;    % Degrees of freedom per node (x,y)
    data.D = D_matrix_Stress(data);
    
    % Setup mesh
    meshSize = meshSizes(i);
    numElemsX = ceil(problemWidth/meshSize);
    numElemsY = ceil(problemHeight/meshSize);
    
    % Adjust mesh size to fit domain exactly
    actualMeshSizeX = problemWidth/numElemsX;
    actualMeshSizeY = problemHeight/numElemsY;
    
    % Calculate number of nodes
    numNodesX = numElemsX + 1;
    numNodesY = numElemsY + 1;
    totalNodes = numNodesX * numNodesY;
    
    % Generate node coordinates
    nodeCoords = zeros(totalNodes, 2);
    nodeCounter = 1;
    
    for j = 1:numNodesY
        for k = 1:numNodesX
            x = (k-1) * actualMeshSizeX;
            y = (j-1) * actualMeshSizeY;
            nodeCoords(nodeCounter, :) = [x, y];
            nodeCounter = nodeCounter + 1;
        end
    end
    
    % Define element connectivity (for 4-noded elements)
    elemConnectivity = zeros(numElemsX * numElemsY, 4);
    elemCounter = 1;
    
    for j = 1:numElemsY
        for k = 1:numElemsX
            n1 = (j-1)*numNodesX + k;           % Bottom left
            n2 = (j-1)*numNodesX + (k+1);       % Bottom right
            n3 = j*numNodesX + (k+1);           % Top right
            n4 = j*numNodesX + k;               % Top left
            
            elemConnectivity(elemCounter, :) = [n1, n2, n3, n4];
            elemCounter = elemCounter + 1;
        end
    end
    
    % Create figure to visualize the mesh
    figure;
    hold on;
    title(sprintf('Mesh with element size = %f', meshSize));
    xlabel('X-coordinate');
    ylabel('Y-coordinate');
    grid on;
    
    % Plot each element with coordinate and node labels
    for elemIdx = 1:size(elemConnectivity, 1)
        elementNodes = elemConnectivity(elemIdx, :);
        nodeCoordinates = nodeCoords(elementNodes, :);
        
        % Close the element by repeating the first node
        nodeCoordinates = [nodeCoordinates; nodeCoordinates(1, :)];
        
        % Plot element with specific color
        plot(nodeCoordinates(:, 1), nodeCoordinates(:, 2), 'b-', 'LineWidth', 1);
        
        % Label the nodes if not too many elements
        if meshSize >= 0.5
            for j = 1:length(elementNodes)
                text(nodeCoords(elementNodes(j), 1), nodeCoords(elementNodes(j), 2), ...
                    num2str(elementNodes(j)), 'FontSize', 8, 'Color', 'r');
            end
        end
    end
    
    % Save the mesh figure
    saveas(gcf, sprintf('mesh_size_%f.png', meshSize));
    
    % Store mesh in data structure
    data.GlobalNodes = nodeCoords;
    data.Connectivity = elemConnectivity;
    data.numElements = size(elemConnectivity, 1);
    data.nTotalDOF = size(nodeCoords, 1) * data.nDOFPNode;
    data.type_of_mesh = "Undistorted";
    data.updateNodalLocations = nodeCoords;
    
    % Set boundary conditions
    % Find nodes on the bottom edge (y = 0)
    bottomNodes = find(abs(nodeCoords(:, 2)) < 1e-6);
    
    % Find nodes on the right edge (x = problemWidth)
    rightNodes = find(abs(nodeCoords(:, 1) - problemWidth) < 1e-6);
    
    % Find top-right corner node (for concentrated force)
    topRightCorner = find(abs(nodeCoords(:, 1) - problemWidth) < 1e-6 & ...
                          abs(nodeCoords(:, 2) - problemHeight) < 1e-6);
    
    % Output top-right corner node info
    fprintf('Force applied at node %d (x=%.2f, y=%.2f)\n', ...
            topRightCorner, nodeCoords(topRightCorner, 1), nodeCoords(topRightCorner, 2));
    
    % Initialize BC arrays
    numBCs = length(bottomNodes) + length(rightNodes);
    data.BCnodes = zeros(numBCs, 1);
    data.BCDof = zeros(numBCs, 1);
    data.BC = zeros(numBCs, 1);
    
    % Set boundary conditions
    bcCounter = 1;
    
    % Bottom edge: restrict y-direction
    for j = 1:length(bottomNodes)
        data.BCnodes(bcCounter) = bottomNodes(j);
        data.BCDof(bcCounter) = 2;  % y-direction constraint
        data.BC(bcCounter) = 0;     % zero displacement
        bcCounter = bcCounter + 1;
    end
    
    % Right edge: restrict x-direction
    for j = 1:length(rightNodes)
        data.BCnodes(bcCounter) = rightNodes(j);
        data.BCDof(bcCounter) = 1;  % x-direction constraint
        data.BC(bcCounter) = 0;     % zero displacement
        bcCounter = bcCounter + 1;
    end
    
    % Set no body forces
    data.f = [0, 0, 0, 0];
    
    % Create elements structure
    elements_structure = create_elements_structure(data);
    
    % Apply concentrated force at top-right corner (F = 20)
    forceValue = 20;
    
    % Create force vector
    F_point = zeros(data.nTotalDOF, 1);
    
    % Apply concentrated force at top-right corner (in y-direction)
    F_point_index = (topRightCorner-1)*data.nDOFPNode + 2; % y-direction DOF
    F_point(F_point_index) = forceValue;
    
    % Create global stiffness matrix and force vector
    Global_elements_structure = Global_Stifness_Matrix(@convert_to_index, elements_structure, data);
    
    % Add concentrated force to global force vector
    Global_elements_structure.F = Global_elements_structure.F + F_point;
    
    % Apply boundary conditions
    [K, F] = assemblyBC(data, Global_elements_structure);
    
    % Solve the system (K*Q = F)
    Q = K \ F;
    Global_elements_structure.Q = Q;
    
    % Calculate displacements and stresses
    Displacementandstraincalc_results = Displacementandstraincalc(Global_elements_structure, data);
    
    % Plot sigma_yy field
    figure;
    hold on;
    title(sprintf('\\sigma_{yy} field with element size = %f', meshSize));
    xlabel('X-coordinate');
    ylabel('Y-coordinate');
    colormap(jet);
    
    % Find min/max stress values
    minStress = Inf;
    maxStress = -Inf;
    
    % Loop through elements to extract and plot stress data
    for elemIdx = 1:data.numElements
        % Get element coordinates
        x = Displacementandstraincalc_results(elemIdx).x;
        y = Displacementandstraincalc_results(elemIdx).y;
        
        % Get stress data - use component 2 for sigma_yy
        stressComponent_values = Displacementandstraincalc_results(elemIdx).Stress(:,:,2);
        
        % Update min/max stress values
        minStress = min(minStress, min(stressComponent_values(:)));
        maxStress = max(maxStress, max(stressComponent_values(:)));
        
        % Create filled contour for this element
        contourf(x, y, stressComponent_values, 20, 'LineStyle', 'none');
    end
    
    % Add colorbar
    colorbar;
    caxis([minStress maxStress]);
    
    % Add element outlines for clarity
    for elemIdx = 1:data.numElements
        elementNodes = elemConnectivity(elemIdx, :);
        nodeCoordinates = nodeCoords(elementNodes, :);
        
        % Close the element by repeating the first node
        nodeCoordinates = [nodeCoordinates; nodeCoordinates(1, :)];
        
        % Plot element outline
        plot(nodeCoordinates(:, 1), nodeCoordinates(:, 2), 'k-', 'LineWidth', 0.5);
    end
    
    % Add text with max/min values
    text(0.02*problemWidth, 0.02*problemHeight, sprintf('Min: %.2f', minStress), ...
         'Units', 'data', 'FontSize', 10, 'BackgroundColor', 'white');
    text(0.02*problemWidth, 0.06*problemHeight, sprintf('Max: %.2f', maxStress), ...
         'Units', 'data', 'FontSize', 10, 'BackgroundColor', 'white');
    
    % Make axis equal for proper visualization
    axis equal;
    axis([0, problemWidth, 0, problemHeight]);
    
    % Save the stress plot
    saveas(gcf, sprintf('sigma_yy_size_%f.png', meshSize));
    
    % Store results for later analysis
    modelResults{i}.data = data;
    modelResults{i}.elements_structure = elements_structure;
    modelResults{i}.Global_elements_structure = Global_elements_structure;
    modelResults{i}.K = K;
    modelResults{i}.F = F;
    modelResults{i}.Q = Q;
    modelResults{i}.Displacementandstraincalc = Displacementandstraincalc_results;
    modelResults{i}.meshSize = meshSize;
    modelResults{i}.maxStress = maxStress;
    modelResults{i}.minStress = minStress;
end

% Save all results
save('singularity_results.mat', 'modelResults', 'meshSizes');

fprintf('\nAnalysis complete. Results saved to singularity_results.mat\n');

% Plot maximum stress values vs mesh size
figure;
maxStressValues = zeros(length(meshSizes), 1);
for i = 1:length(meshSizes)
    maxStressValues(i) = modelResults{i}.maxStress;
end

loglog(meshSizes, maxStressValues, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Element Size');
ylabel('Maximum \sigma_{yy}');
title('Maximum \sigma_{yy} vs Element Size');
grid on;
saveas(gcf, 'max_stress_vs_mesh_size.png');

% Add this after all the individual stress plots but before saving the results

% Create a new figure for plotting sigma_yy vs y at x=2
figure;
hold on;
title('\sigma_{yy} along x=2 for different mesh sizes');
xlabel('y-coordinate');
ylabel('\sigma_{yy}');
grid on;

% Line styles and colors for different mesh sizes
lineStyles = {'-o', '-s', '-^', '-d'};
lineColors = {'b', 'r', 'g', 'm'};
legendEntries = cell(length(meshSizes), 1);

% Loop through each mesh size
for i = 1:length(meshSizes)
    % Get the mesh size for legend
    meshSize = meshSizes(i);
    legendEntries{i} = sprintf('Element size = %.3f', meshSize);
    
    % Extract data for this mesh
    Displacementandstraincalc_results = modelResults{i}.Displacementandstraincalc;
    nodeCoords = modelResults{i}.data.GlobalNodes;
    
    % Find nodes along x=2 (with small tolerance for numerical precision)
    xLine = 2.0;
    x2Nodes = find(abs(nodeCoords(:,1) - xLine) < 1e-6);
    
    % Get y-coordinates and sort them
    y2Coords = nodeCoords(x2Nodes, 2);
    [y2Coords, sortIdx] = sort(y2Coords);
    x2Nodes = x2Nodes(sortIdx);
    
    % Initialize arrays for storing stress values
    sigma_yy_values = zeros(size(y2Coords));
    
    % Get right edge elements (elements that contain nodes on x=2)
    elements = modelResults{i}.data.Connectivity;
    numElements = modelResults{i}.data.numElements;
    
    % Loop through nodes at x=2 to extract stress values
    for j = 1:length(x2Nodes)
        node = x2Nodes(j);
        yVal = y2Coords(j);
        
        % Find which element contains this node
        containingElements = [];
        for elemIdx = 1:numElements
            if any(elements(elemIdx,:) == node)
                containingElements = [containingElements, elemIdx];
            end
        end
        
        % Average the stresses from elements that contain this node
        if ~isempty(containingElements)
            stressSum = 0;
            for elemIdx = containingElements
                % Get the element's stress field for sigma_yy
                stressField = Displacementandstraincalc_results(elemIdx).Stress(:,:,2);
                
                % Get coordinates of the element
                elemX = Displacementandstraincalc_results(elemIdx).x;
                elemY = Displacementandstraincalc_results(elemIdx).y;
                
                % Find the closest point in the element's stress field to this node
                [~, minDistIdx] = min(sqrt((elemX(:) - xLine).^2 + (elemY(:) - yVal).^2));
                [row, col] = ind2sub(size(elemX), minDistIdx);
                
                % Add the stress value at this point
                stressSum = stressSum + stressField(row, col);
            end
            
            % Calculate average stress
            sigma_yy_values(j) = stressSum / length(containingElements);
        end
    end
    
    % Plot sigma_yy vs y for this mesh size
    plot(y2Coords, sigma_yy_values, lineStyles{i}, 'Color', lineColors{i}, 'LineWidth', 1.5, 'MarkerSize', 6);
end

% Add legend
legend(legendEntries, 'Location', 'best');

% Add additional features
grid on;
title('\sigma_{yy} along x=2 for different mesh sizes');
xlabel('y-coordinate');
ylabel('\sigma_{yy}');

% Save the figure
saveas(gcf, 'sigma_yy_vs_y_at_x2.png');