function [Strain] = strain_cal(B, q)
% strain_cal - Calculates strain using B matrix and displacement vector
%
% PSEUDOCODE:
% 1. Multiply the B matrix by the displacement vector q
% 2. Return the resulting strain vector
%
% Input:
%   B - B matrix relating displacement derivatives to strain
%   q - Displacement vector (nodal displacements)
%
% Output:
%   Strain - Calculated strain vector

Strain = B * q;                       % Calculate strain using the relationship: ε = B·q
end