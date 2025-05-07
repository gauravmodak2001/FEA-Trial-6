function updatedNodeCoordinates = updateNodalLocations(data)
% updateNodalLocations - Updates node coordinates for distorted/undistorted mesh
%
% PSEUDOCODE:
% 1. Extract the original node coordinates from input data structure
% 2. Calculate number of nodes to be distorted based on corner nodes
% 3. Determine mesh size and type (distorted or undistorted)
% 4. Initialize the output coordinates with the original values
% 5. If mesh type is "distorted":
%    a. Loop through the nodes that need to be distorted
%    b. For nodes with odd indices, move them upward (increase y-coordinate)
%    c. For nodes with even indices, move them downward (decrease y-coordinate)
% 6. If mesh type is "Undistorted", keep original coordinates
% 7. Return the updated node coordinates
%
% Input:
%   data - Structure containing mesh information including node coordinates
%
% Output:
%   updatedNodeCoordinates - Updated node coordinates after distortion

nodeCoordinates = data.GlobalNodes;    % Extract original node coordinates from input data
distorted_nodes = data.Total_corner_nodes - (data.top_corner_node + data.bottom_corner_node); % Calculate number of nodes to be distorted
h = data.mesh_size;                   % Extract mesh size parameter
type = (data.type_of_mesh);           % Extract mesh type information
updatedNodeCoordinates = nodeCoordinates; % Initialize output with original coordinates

if type == "distorted"                % Check if mesh type is distorted
    for i = data.bottom_corner_node+1:(distorted_nodes+data.bottom_corner_node) % Loop through nodes that need distortion
        if mod(i, 2) == 1             % Check if node index is odd
            updatedNodeCoordinates(i, 2) = nodeCoordinates(i, 2) + data.d * h / 2; % Move odd-indexed nodes upward
        else                          % If node index is even
            updatedNodeCoordinates(i, 2) = nodeCoordinates(i, 2) - data.d * h / 2; % Move even-indexed nodes downward
        end
    end
elseif type == "Undistorted"          % Check if mesh type is undistorted
    updatedNodeCoordinates = data.GlobalNodes; % Keep original coordinates for undistorted mesh
end
end