function[]= xy_coordinate_plot(elements_structure)
% xy_coordinate_plot - Visualizes shape functions over x-y coordinates
%
% PSEUDOCODE:
% 1. Extract shape functions and coordinate data from the input structure
% 2. Determine dimensions of the shape function matrix
% 3. Set up subplot layout based on number of shape functions
% 4. Loop through shape functions (odd indices only to skip zero values)
%    a. Create a subplot for each shape function
%    b. Extract the shape function values at the current index
%    c. Create a filled contour plot of the shape function
%    d. Add labels, colorbar, and title
%    e. Adjust axis settings
%
% Input:
%   elements_structure - Structure containing shape functions and coordinates
%
% Output:
%   Visualization of shape functions as contour plots (no return value)

%% Plotting x and y vs the Nodes %% 
N = elements_structure.shapeFunctions; % Extract the shape functions matrix from input structure
x = elements_structure.x;              % Extract the x-coordinate matrix from input structure
y = elements_structure.y;              % Extract the y-coordinate matrix from input structure

% Check for the dimensions of matrix N
[n, m, p] = size(N);                   % Get dimensions: n rows, m columns, p is depth (number of shape functions)

% Divide the number of plots equally for subplots
l = p/2;                               % Calculate number of columns for subplot grid (half of total shape functions)

% Loop through the plots with a step size of 2 to avoid zero values
for plotId = 1:2:p                     % Loop with step size 2, starting from 1 (1, 3, 5, ...) to skip zero-valued functions
    
    % Create subplot
    subplot(2, l, plotId);             % Create a subplot in a 2×l grid at position plotId
    
    % Access the 3rd dimension of matrix to plot
    z = N(:, :, plotId);               % Extract the specific shape function at index plotId from the 3D matrix
    
    % Plotting commands and options
    contourf(x, y, z);                 % Create filled contour plot with x and y coordinates and z values
    xlabel(['x']);                     % Add label for x-axis
    ylabel(['y']);                     % Add label for y-axis
    colorbar;                          % Add colorbar to show mapping of colors to function values
    title(strcat('N', num2str(plotId))); % Set title as 'N1', 'N3', etc. based on current plotId
    axis auto;                         % Automatically set axis limits based on data range
end
end