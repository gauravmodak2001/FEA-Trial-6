function plotMesh(model)
    % Function to plot the mesh (nodes and elements)
    % Input:
    %   model - structure containing mesh data
    
    % Get node coordinates
    nodes = model.GlobalNodes;
    
    % Get element connectivity
    connectivity = model.Connectivity;
    
    % Create figure
    hold on;
    title('Element Connectivity Visualization');
    xlabel('X-coordinate');
    ylabel('Y-coordinate');
    grid on;
    
    % Define colors for elements (cycling through a set of colors)
    colorSet = ['r', 'g', 'b', 'c', 'm', 'y'];
    
    % Plot each element
    for i = 1:model.numElements
        % Get the 4 corner nodes (ignoring duplicates in 8-node format)
        elementNodes = connectivity(i, 1:4);
        nodeCoords = nodes(elementNodes, :);
        
        % Close the element by repeating the first node
        nodeCoords = [nodeCoords; nodeCoords(1, :)];
        
        % Plot element with specific color
        colorIdx = mod(i-1, length(colorSet)) + 1;
        plot(nodeCoords(:, 1), nodeCoords(:, 2), colorSet(colorIdx), 'LineWidth', 1);
        
        % Calculate element center for element number
        centerX = mean(nodeCoords(1:end-1, 1));
        centerY = mean(nodeCoords(1:end-1, 2));
        text(centerX, centerY, sprintf('%d', i), 'FontSize', 8, 'Color', 'k', 'HorizontalAlignment', 'center');
    end
    
    % Plot and label all nodes
    scatter(nodes(:, 1), nodes(:, 2), 30, 'k', 'filled');
    for j = 1:size(nodes, 1)
        text(nodes(j, 1), nodes(j, 2), sprintf(' %d', j), 'FontSize', 8, 'Color', 'b');
    end
    
    % Add additional annotations
    title(sprintf('Mesh with %d elements and %d nodes', model.numElements, size(nodes, 1)));
    
    % Set axis limits with some padding
    maxWidth = max(nodes(:, 1));
    maxHeight = max(nodes(:, 2));
    axis([-0.1, maxWidth*1.1, -0.1, maxHeight*1.1]);
    
    % Make axis equal for proper visualization
    axis equal;
    hold off;
end