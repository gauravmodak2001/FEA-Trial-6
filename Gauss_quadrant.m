function [Integration_Points, Weights] = Gauss_quadrant(num_gauss_points)
% PSEUDOCODE:
% 1. Define tables of Gauss points and weights for different orders
% 2. Based on requested number of points, select appropriate points and weights
% 3. Return selected integration points and weights

% Function to calculate integration points using Gaussian quadrature
% Table of Gauss points from 1-point to 6-point rules
XiEta_coor = [0, 0, 0, 0, 0, 0;                                       % 1-point rule
            -0.57735, 0.57735, 0, 0, 0, 0;                           % 2-point rule
             0, 0.77460, -0.77460, 0, 0, 0;                          % 3-point rule
             0.861136, 0.339981, -0.339981, -0.861136, 0, 0;         % 4-point rule
             0.906179, 0.538469, 0, -0.538469, -0.906179, 0;         % 5-point rule
             0.93246, 0.661209, 0.23861, -0.23861, -0.661209, -0.93246]; % 6-point rule

% Table of weights corresponding to Gauss points
W = [2, 0, 0, 0, 0, 0;                                                % 1-point rule
    1.0, 1.0, 0, 0, 0, 0;                                            % 2-point rule
    0.88889, 0.55556, 0.55556, 0, 0, 0;                              % 3-point rule
    0.347854, 0.652145, 0.652145, 0.347854, 0, 0;                    % 4-point rule
    0.236926, 0.478628, 0.56888, 0.478628, 0.236926, 0;              % 5-point rule
    0.1713244, 0.36076, 0.467913, 0.467913, 0.36076, 0.1713244];     % 6-point rule

% Select appropriate points and weights based on requested number
if num_gauss_points == 1
   % Single integration point: ∫f(x)dx ≈ w₁f(x₁)
   Integration_Points = [XiEta_coor(1,1)];
   Weights = [W(1,1)];
elseif num_gauss_points == 2
   % Two-point rule: ∫f(x)dx ≈ w₁f(x₁) + w₂f(x₂)
   Integration_Points = [XiEta_coor(2,1), XiEta_coor(2,2)];
   Weights = [W(2,1), W(2,2)];
elseif num_gauss_points == 3
   % Three-point rule: ∫f(x)dx ≈ w₁f(x₁) + w₂f(x₂) + w₃f(x₃)
   Integration_Points = [XiEta_coor(3,1), XiEta_coor(3,2), XiEta_coor(3,3)];
   Weights = [W(3,1), W(3,2), W(3,3)];
elseif num_gauss_points == 4
   % Four-point rule
   Integration_Points = [XiEta_coor(4,1), XiEta_coor(4,2), XiEta_coor(4,3), XiEta_coor(4,4)];
   Weights = [W(4,1), W(4,2), W(4,3), W(4,4)];
elseif num_gauss_points == 5
   % Five-point rule
   Integration_Points = [XiEta_coor(5,1), XiEta_coor(5,2), XiEta_coor(5,3), XiEta_coor(5,4), XiEta_coor(5,5)];
   Weights = [W(5,1), W(5,2), W(5,3), W(5,4), W(5,5)];
elseif num_gauss_points == 6
   % Six-point rule
   Integration_Points = [XiEta_coor(6,1), XiEta_coor(6,2), XiEta_coor(6,3), XiEta_coor(6,4), XiEta_coor(6,5), XiEta_coor(6,6)];
   Weights = [W(6,1), W(6,2), W(6,3), W(6,4), W(6,5), W(6,6)];
else
   error('Unsupported number of Gauss points. Use 1, 2 or 3');
end
end