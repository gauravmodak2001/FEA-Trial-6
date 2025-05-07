function [] = Displacement_plot(Displacementandstraincalc, data)
% PSEUDOCODE:
% 1. For each element in the model:
%    a. Extract displacement data and coordinates
%    b. For each displacement component (u, v):
%       i.   Create a subplot
%       ii.  Extract component values
%       iii. Generate contour plot
%       iv.  Add labels and formatting
% 2. Add a global title for the figure

%% Plotting Displacement vs Xi and Eta %%
% Loop through the elements
for i = 1:data.numElements
   % Access displacement data for the i-th element
   u = Displacementandstraincalc(i).U;          % Get displacement field
   x_values = Displacementandstraincalc(i).x;   % Get x-coordinates
   y_values = Displacementandstraincalc(i).y;   % Get y-coordinates
   
   % Check the size of displacement data
   [~, ~, p] = size(u);                         % p is the number of components (usually 2 for 2D)
   
   % Loop through the displacement components
   for plotId = 1:p
       % Create subplot
       subplot(1, p, plotId);                   % Create 1×p grid of subplots
       
       % Extract values for plotting
       w = u(:, :, plotId);                     % Extract displacement component (u₁ or u₂)
       
       % Plot the data
       contourf(x_values, y_values, w);         % Create filled contour plot
       
       % Adjust mesh grid for smoother appearance
       shading interp;                          % Use interpolated shading for smooth plot
       colormap jet;                            % Use jet colormap
       
       % Add labels and title
       xlabel('x');                             % x-axis label
       ylabel('y');                             % y-axis label
       zlabel(['u', num2str(plotId)]);          % z-axis label (u₁ or u₂)
       title(['Displacement Component ', num2str(plotId)]); % Component title
       
       % Hold on for subsequent plots
       hold on;
   end
end

% Set global title for entire figure
sgtitle('Displacements');                        % Add super title for the figure
end