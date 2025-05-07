function [index] = convert_to_index(Node_num, Dir, nDOFPNODE)

index = nDOFPNODE * (Node_num - 1) + Dir;  % Calculate global DOF index
end