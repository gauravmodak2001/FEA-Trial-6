function elements_structure = createElementsStructure(model)
    % Create elements structure for FEA analysis
    % Adapted for 4-node elements
    
    connectivityMatrix = model.Connectivity;
    allNodalLocations = model.GlobalNodes;
    D = model.D;
    nDOFPNode = model.nDOFPNode;
    
    numElements = size(connectivityMatrix, 1);
    
    % Initialize output structure
    elements_structure(numElements).shapeFunctions = [];
    elements_structure(numElements).DiffShapefunctions = [];
    elements_structure(numElements).Stiffnessmatrix = [];
    elements_structure(numElements).SurfaceForce = [];
    elements_structure(numElements).x = [];
    elements_structure(numElements).y = [];
    
    % Setup Gaussian quadrature points and weights
    Gauss_in_Xi = 2;  % Number of Gauss points in ξ direction
    Gauss_in_Eta = 2; % Number of Gauss points in η direction
    [Integration_points, Weights] = Gauss_quadrant(2); % Using 2-point rule for 4-node elements
    
    % Loop through each element
    for i = 1:numElements
        % Get element data
        nodeIndices = connectivityMatrix(i, :);
        X = allNodalLocations(nodeIndices, :);
        X = reshape(X.', [], 1);
        
        % Calculate element stiffness matrix
        K_temp = Calculate_K_matrix(4, X, D, Gauss_in_Xi, Gauss_in_Eta, Integration_points, Weights, model);
        
        % Initialize surface force vector
        SurfaceForce = zeros(nDOFPNode * 4, 1); % For 4-node elements
        
        % Calculate shape functions and coordinates at integration points
        n = 10; % Grid size for visualization
        x = zeros(n+1, n+1);
        y = zeros(n+1, n+1);
        
        for j = 1:n+1
            Xi_j = -1 + 2*(j-1)/n;
            for k = 1:n+1
                Eta_k = -1 + 2*(k-1)/n;
                
                % Calculate shape functions
                shapeFunctions = Shape_function([Xi_j, Eta_k], 4);
                DiffShapefunctions = Diff_Shape_functions([Xi_j, Eta_k], 4);
                
                % Calculate physical coordinates
                temp_C = calculate_coordinates(@Shape_function, [Xi_j, Eta_k], X, 4);
                x(j, k) = temp_C(1, 1);
                y(j, k) = temp_C(2, 1);
            end
        end
        
        % Store element properties in output structure
        elements_structure(i).SurfaceForce = SurfaceForce;
        elements_structure(i).shapeFunctions = shapeFunctions;
        elements_structure(i).DiffShapefunctions = DiffShapefunctions;
        elements_structure(i).Stiffnessmatrix = K_temp;
        elements_structure(i).x = x;
        elements_structure(i).y = y;
    end
end