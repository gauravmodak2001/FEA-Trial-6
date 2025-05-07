function plotStressField(results, model, stressComponent)
    % Function to plot stress field
    % Inputs:
    %   results - structure with stress/strain results from Displacementandstraincalc
    %   model - structure containing mesh data
    %   stressComponent - index of stress component to plot (1=σxx, 2=σyy, 3=τxy)
    
    % Create a new figure
    hold on;
    
    % Define colormap
    colormap(jet);
    
    % Initialize arrays to store stress values for contour plot
    nodes = model.GlobalNodes;
    numNodes = size(nodes, 1);
    xCoords = nodes(:,1);
    yCoords = nodes(:,2);
    stressValues = zeros(numNodes, 1);
    
    % Store min/max stress for colorbar
    minStress = Inf;
    maxStress = -Inf;
    
    % Create triangulation for elements
    elements = model.Connectivity;
    numElements = model.numElements;
    
    % Loop through elements to extract and plot stress data
    for i = 1:numElements
        % Get stress data for current element
        elementStress = results(i).Stress;
        
        % Extract the component we want
        stressComponent_values = elementStress(:,:,stressComponent);
        
        % Update min/max stress
        minStress = min(minStress, min(stressComponent_values(:)));
        maxStress = max(maxStress, max(stressComponent_values(:)));
        
        % Get element coordinates
        x = results(i).x;
        y = results(i).y;
        
        % Create filled contour for this element
        contourf(x, y, stressComponent_values, 20, 'LineStyle', 'none');
    end
    
    % Add colorbar
    colorbar;
    caxis([minStress maxStress]);
    
    % Add element outlines for clarity
    for i = 1:numElements
        % Get the 4 corner nodes (for 4-node elements)
        % Note: If using 8-node connectivity matrix, only use corners
        elementNodes = elements(i, 1:4);
        nodeCoords = nodes(elementNodes, :);
        
        % Close the element by repeating the first node
        nodeCoords = [nodeCoords; nodeCoords(1, :)];
        
        % Plot element outline
        plot(nodeCoords(:, 1), nodeCoords(:, 2), 'k-', 'LineWidth', 0.5);
    end
    
    % Add title and labels
    if stressComponent == 1
        title('\sigma_{xx} Stress Distribution');
    elseif stressComponent == 2
        title('\sigma_{yy} Stress Distribution');
    elseif stressComponent == 3
        title('\tau_{xy} Stress Distribution');
    end
    xlabel('X-coordinate');
    ylabel('Y-coordinate');
    
    % Set axis limits
    maxWidth = max(nodes(:, 1));
    maxHeight = max(nodes(:, 2));
    axis([0, maxWidth, 0, maxHeight]);
    
    % Make axis equal for proper visualization
    axis equal;
    grid on;
    
    % Add text with max/min values
    text(0.02*maxWidth, 0.02*maxHeight, sprintf('Min: %.2f', minStress), ...
         'Units', 'data', 'FontSize', 10, 'BackgroundColor', 'white');
    text(0.02*maxWidth, 0.06*maxHeight, sprintf('Max: %.2f', maxStress), ...
         'Units', 'data', 'FontSize', 10, 'BackgroundColor', 'white');
    
    hold off;
end