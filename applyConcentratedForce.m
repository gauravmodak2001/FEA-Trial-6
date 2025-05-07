function model = applyConcentratedForce(model, forceValue)
    % Apply concentrated force at top-right corner
    
    % Find the top-right corner node (max x, max y)
    nodes = model.GlobalNodes;
    [maxX, ~] = max(nodes(:,1));
    [maxY, ~] = max(nodes(:,2));
    
    % Find node closest to top-right corner
    distances = sqrt((nodes(:,1) - maxX).^2 + (nodes(:,2) - maxY).^2);
    [~, cornerNodeIndex] = min(distances);
    
    fprintf('Applying force at node %d (x=%.2f, y=%.2f)\n', ...
            cornerNodeIndex, nodes(cornerNodeIndex,1), nodes(cornerNodeIndex,2));
    
    % Initialize force fields if they don't exist
    if ~isfield(model, 'forceNodes')
        model.forceNodes = [];
    end
    if ~isfield(model, 'forceDof')
        model.forceDof = [];
    end
    if ~isfield(model, 'forceValues')
        model.forceValues = [];
    end
    
    % Add this force to the model
    model.forceNodes = [model.forceNodes; cornerNodeIndex];
    model.forceDof = [model.forceDof; 2]; % y-direction
    model.forceValues = [model.forceValues; forceValue];
end