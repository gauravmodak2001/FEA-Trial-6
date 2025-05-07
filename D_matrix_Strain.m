function D_pl_strain = D_matrix_Strain(data)
% PSEUDOCODE:
% 1. Extract material properties from data structure
% 2. Calculate D matrix for plane strain condition
% 3. Return D matrix

% Calculate D matrix for Plane Strain condition
E = data.E;                           % Young's modulus
v = data.v;                           % Poisson's ratio

% Calculate constitutive matrix D for plane strain condition
% Formula: D = E/((1+v)(1-2v)) * [1-v    v      0    ]
%                                [v      1-v    0    ]
%                                [0      0    (1-2v)/2]
D_pl_strain = (E / ((1 + v) * (1 - 2 * v))) * [1 - v, v, 0; 
                                             v, 1 - v, 0; 
                                             0, 0, 0.5 * (1 - 2 * v)];
end