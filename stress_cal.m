function [stress] = stress_cal(strain_cal, B, q, D_pl_stress)
% stress_cal - Calculates stress from strain using constitutive relationship
%
% PSEUDOCODE:
% 1. Calculate strain using the strain_cal function
% 2. Calculate stress by multiplying the D matrix with the strain vector
% 3. Return the stress vector
%
% Input:
%   strain_cal - Function handle to calculate strain
%   B - B matrix relating displacement derivatives to strain
%   q - Displacement vector
%   D_pl_stress - Material property matrix (constitutive matrix)
%
% Output:
%   stress - Calculated stress vector

strain = strain_cal(B, q);            % Calculate strain using the provided strain_cal function
stress = D_pl_stress * strain;        % Calculate stress using constitutive relationship σ = D·ε
end