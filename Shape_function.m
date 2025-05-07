function [N] = Shape_function(Xi_coor, s)
% Shape_function - Calculates shape functions in parametric coordinates
%
% PSEUDOCODE:
% 1. Extract Xi and Eta coordinates from input parameter
% 2. Based on number of nodes parameter s:
%    a. If s=8, calculate shape functions for 8-node quadrilateral element
%    b. If s=4, calculate shape functions for 4-node quadrilateral element
%    c. Otherwise, throw an error for unsupported node configurations
% 3. Arrange shape functions into a matrix format for vector calculations
% 4. Return the shape function matrix
%
% Input:
%   Xi_coor - Matrix containing parametric coordinates [Xi, Eta]
%   s - Number of nodes in the element (4 or 8)
%
% Output:
%   N - Matrix of shape functions arranged for displacement interpolation

% Extract coordinates for clarity
Xi = Xi_coor(:,1);                    % Extract Xi coordinate from input
Eta = Xi_coor(:,2);                   % Extract Eta coordinate from input

% Calculate shape functions based on the number of nodes
if s == 4                        % Check if using 4-node element
    % Shape functions for 4-node configuration (bilinear element)
    N1 = (Xi-1)*(Eta-1)*0.25;              % Corner node 1 shape function
    N2 = (Xi+1)*(Eta-1)*(-0.25);           % Corner node 2 shape function
    N3 = (Xi+1)*(Eta+1)*0.25;              % Corner node 3 shape function
    N4 = (Xi-1)*(Eta+1)*(-0.25);           % Corner node 4 shape function
    
    % Arrange shape functions into a matrix for vector calculations
    N = [N1,0,N2,0,N3,0,N4,0;              % Row for x-direction (u)
         0,N1,0,N2,0,N3,0,N4];             % Row for y-direction (v)
         
elseif s == 8                            % Check if using 8-node element
    % Shape functions for 8-node configuration (quadratic serendipity element)
    N1 = (1-Xi)*(1-Eta)*(1+Xi+Eta)*(-0.25);  % Corner node 1 shape function
    N2 = (1-Xi^2)*(1-Eta)*0.5;               % Mid-side node 2 shape function
    N3 = (1+Xi)*(1-Eta)*(1-Xi+Eta)*(-0.25);  % Corner node 3 shape function
    N4 = (1+Xi)*(1-Eta^2)*0.5;               % Mid-side node 4 shape function
    N5 = (1+Xi)*(1+Eta)*(1-Xi-Eta)*(-0.25);  % Corner node 5 shape function
    N6 = (1-Xi^2)*(1+Eta)*0.5;               % Mid-side node 6 shape function
    N7 = (1-Xi)*(1+Eta)*(1+Xi-Eta)*(-0.25);  % Corner node 7 shape function
    N8 = (1-Xi)*(1-Eta^2)*0.5;               % Mid-side node 8 shape function
    
    % Arrange shape functions into a matrix for vector calculations
    % Each row corresponds to a displacement direction (u,v)
    % Each column pair corresponds to a node's two DOFs
    N = [N1,0,N2,0,N3,0,N4,0,N5,0,N6,0,N7,0,N8,0;  % Row for x-direction (u)
         0,N1,0,N2,0,N3,0,N4,0,N5,0,N6,0,N7,0,N8]; % Row for y-direction (v)
         
         
else
    % Handle unsupported node configurations
    error('Unsupported number of nodes. Use either 4 or 8 nodes.This function is not supported for 2 and 3 noded element');  % Throw error for invalid input
end
end