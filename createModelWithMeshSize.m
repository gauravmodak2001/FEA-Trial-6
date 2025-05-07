function model = createModelWithMeshSize(meshSize, width, height, E, v, thickness)
    % Create a model with specified mesh size and problem dimensions
    
    % Initialize model data structure
    model = struct();
    model.E = E;
    model.v = v;
    model.thikness = thickness;  % Note: keep the original spelling "thikness" for compatibility
    model.nDOFPNode = 2;    % Degrees of freedom per node (x,y)
    
    % Define plane stress condition
    model.D = D_matrix_Stress(model);
    
    % Calculate number of elements in each direction
    numElemsX = ceil(width/meshSize);
    numElemsY = ceil(height/meshSize);
    
    % Adjust mesh size to fit domain exactly
    actualMeshSizeX = width/numElemsX;
    actualMeshSizeY = height/numElemsY;
    
    % Calculate number of nodes
    numNodesX = numElemsX + 1;
    numNodesY = numElemsY + 1;
    totalNodes = numNodesX * numNodesY;
    
    % Generate node coordinates
    nodeCoords = zeros(totalNodes, 2);
    nodeCounter = 1;
    
    for j = 1:numNodesY
        for i = 1:numNodesX
            x = (i-1) * actualMeshSizeX;
            y = (j-1) * actualMeshSizeY;
            nodeCoords(nodeCounter, :) = [x, y];
            nodeCounter = nodeCounter + 1;
        end
    end
    
    % Define element connectivity (for 4-noded elements)
    elemConnectivity = zeros(numElemsX * numElemsY, 4);
    elemCounter = 1;
    
    for j = 1:numElemsY
        for i = 1:numElemsX
            n1 = (j-1)*numNodesX + i;           % Bottom left
            n2 = (j-1)*numNodesX + (i+1);       % Bottom right
            n3 = j*numNodesX + (i+1);           % Top right
            n4 = j*numNodesX + i;               % Top left
            
            elemConnectivity(elemCounter, :) = [n1, n2, n3, n4];
            elemCounter = elemCounter + 1;
        end
    end
    
    % Store mesh data in model
    model.GlobalNodes = nodeCoords;
    model.numElements = size(elemConnectivity, 1);
    
    % For 4-node elements, we'll convert to 4-node format in Connectivity
    model.Connectivity = elemConnectivity;
    
    % Set total DOFs
    model.nTotalDOF = size(model.GlobalNodes, 1) * model.nDOFPNode;
    model.type_of_mesh = "Undistorted";
    model.updateNodalLocations = model.GlobalNodes;
    
    % Define boundary conditions
    % Find nodes on the bottom edge (y = 0)
    bottomNodes = find(abs(nodeCoords(:, 2)) < 1e-6);
    
    % Find nodes on the right edge (x = width)
    rightNodes = find(abs(nodeCoords(:, 1) - width) < 1e-6);
    
    % Initialize BC arrays
    totalBCs = length(bottomNodes) + length(rightNodes);
    model.BCnodes = zeros(totalBCs, 1);
    model.BCDof = zeros(totalBCs, 1);
    model.BC = zeros(totalBCs, 1);
    
    % Set boundary conditions:
    % - Bottom edge: restrict y-direction (vertical)
    % - Right edge: restrict x-direction (horizontal)
    bcCounter = 1;
    
    for i = 1:length(bottomNodes)
        model.BCnodes(bcCounter) = bottomNodes(i);
        model.BCDof(bcCounter) = 2;  % y-direction constraint
        model.BC(bcCounter) = 0;     % zero displacement
        bcCounter = bcCounter + 1;
    end
    
    for i = 1:length(rightNodes)
        model.BCnodes(bcCounter) = rightNodes(i);
        model.BCDof(bcCounter) = 1;  % x-direction constraint
        model.BC(bcCounter) = 0;     % zero displacement
        bcCounter = bcCounter + 1;
    end
    
    % Store mesh size in model for reference
    model.meshSize = meshSize;
    model.actualMeshSizeX = actualMeshSizeX;
    model.actualMeshSizeY = actualMeshSizeY;
    model.numElemsX = numElemsX;
    model.numElemsY = numElemsY;
    
    % Set no body forces initially
    model.f = [0, 0, 0, 0];
end