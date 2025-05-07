function plot_Stress_Strain_Contours_Xi_Eta(Displacementandstraincalc, data)
% Plot individual contour plots for stress and strain components
% Initialize figure
figure;

% Plot Stress and Strain for each element
for j = 1:data.numElements
    Strain = Displacementandstraincalc(j).Strain; % Get strain data for current element
    Stress = Displacementandstraincalc(j).Stress; % Get stress data for current element
    y = Displacementandstraincalc(j).y;           % Get y coordinates
    x = Displacementandstraincalc(j).x;           % Get x coordinates
    nRows = size(Strain, 3);                      % Number of components
    
    % Plot Stress (S12 only)
    subplot(1, nRows, 1);
    hold on;                                      % Enable hold on for subsequent plots
    surf(x, y, Stress(:,:,2), EdgeColor="black"); % Create surface plot for shear stress S12
    
    % Adjust mesh grid for smoother appearance
    shading interp;                               % Use interpolated shading for smoother plot
    colormap jet;                                 % Use jet colormap
    title('Stress - S12');                        % Add title
    xlabel('x');                                  % Label x-axis
    ylabel('y');                                  % Label y-axis
    zlabel('Stress');                             % Label z-axis
    colorbar                                      % Add colorbar
end
end