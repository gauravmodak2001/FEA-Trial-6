function Global_elements_structure = Global_Stifness_Matrix(convert_to_index, elements_structure, data)
    % Modified version to handle concentrated forces
    
    connectivityMatrix = data.Connectivity;
    nTotalDof = data.nTotalDOF;
    nDOFPNODE = data.nDOFPNode;
    
    % Adjust for 4-node vs 8-node elements
    if size(connectivityMatrix, 2) == 4
        num_elemental_nodes = 4;
    else
        num_elemental_nodes = 8;
    end
    
    numElements = size(connectivityMatrix, 1);
    
    K = zeros(nTotalDof, nTotalDof);
    F = zeros(nTotalDof, 1);
    
    % Initialize output structure
    Global_elements_structure = struct();
    
    % Assembly loop for stiffness matrix
    for i = 1:numElements
        K_elements = elements_structure(i).Stiffnessmatrix;
        
        for j = 1:num_elemental_nodes
            for k = 1:nDOFPNODE
                elemental_row = convert_to_index(j, k, nDOFPNODE);
                Local_row_Node = connectivityMatrix(i, j);
                Local_Row = convert_to_index(Local_row_Node, k, nDOFPNODE);
                
                for l = 1:num_elemental_nodes
                    for m = 1:nDOFPNODE
                        elemental_col = convert_to_index(l, m, nDOFPNODE);
                        Local_Col_Node = connectivityMatrix(i, l);
                        Local_Col = convert_to_index(Local_Col_Node, m, nDOFPNODE);
                        
                        K(Local_Row, Local_Col) = K(Local_Row, Local_Col) + K_elements(elemental_row, elemental_col);
                    end
                end
            end
        end
    end
    
    % Add concentrated forces if specified
    if isfield(data, 'forceNodes') && ~isempty(data.forceNodes)
        for i = 1:length(data.forceNodes)
            node = data.forceNodes(i);
            dof = data.forceDof(i);
            value = data.forceValues(i);
            
            % Calculate global DOF index
            globalDOF = convert_to_index(node, dof, nDOFPNODE);
            
            % Apply force to global force vector
            F(globalDOF) = F(globalDOF) + value;
        end
    end
    
    % Store assembled matrices in output structure
    Global_elements_structure.K = K;
    Global_elements_structure.F = F;
end