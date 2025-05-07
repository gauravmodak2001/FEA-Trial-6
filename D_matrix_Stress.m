function D_pl_stress = D_matrix_Stress(data)
% PSEUDOCODE:
% 1. Extract material properties from data structure
% 2. Calculate D matrix for plane stress condition
% 3. Return D matrix

% Calculate D matrix for Plane Stress condition
E = data.E;                           % Young's modulus
v = data.v;                           % Poisson's ratio

% Calculate constitutive matrix D for plane stress condition
% Formula: D = E/(1-v²) * [1  v  0]
%                         [v  1  0]
%                         [0  0 (1-v)/2]
D_pl_stress = (E / (1 - v^2)) * [1, v, 0; 
                                v, 1, 0; 
                                0, 0, (1 - v)/2];
end