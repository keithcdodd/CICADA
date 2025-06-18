function IC_assignment_table = assign_ICs_to_networks_clust(networks, IC_masks, network_names)
% assign_ICs_to_networks_clust - Assigns ICs to canonical networks by greedy Dice maximization
%
% This function selects which ICs, when combined, best approximate each network mask.
% For each network, it:
%   1. Computes Dice overlap between each IC and the network.
%   2. Greedily adds ICs that improve Dice with the union of selected ICs.
%   3. Stops when 3 consecutive ICs do not improve the Dice score.
%
% Inputs:
%   networks       - 4D binary matrix (X x Y x Z x N), each volume is a canonical network mask
%   IC_masks       - 4D probabilistic matrix (X x Y x Z x M), each volume is a probabilistic IC map
%   network_names  - (Optional) Cell array of strings for each network (length N)
%
% Output:
%   IC_assignment_table - Table with selected ICs, assigned network index, and network name

    num_networks = size(networks, 4);
    num_ICs = size(IC_masks, 4);

    % Generate default network names if none provided
    if nargin < 3 || isempty(network_names)
        network_names = strcat("Network_", string(1:num_networks));
    elseif length(network_names) ~= num_networks
        error('Length of network_names must match number of networks (4th dim of `networks`).');
    end

    % Step 1: Threshold ICs to create binary masks (prob > 0.999)
    binary_ICs = false(size(IC_masks));
    for i = 1:num_ICs
        binary_ICs(:,:,:,i) = IC_masks(:,:,:,i) > 0.999;
    end

    % Initialize arrays to store final IC assignments
    final_ICs = [];                % List of selected IC indices
    final_network_indices = [];   % Corresponding network indices
    final_network_labels = [];    % Corresponding network names

    % Step 2: Process each network independently
    for net = 1:num_networks
        net_mask = logical(networks(:,:,:,net));  % Binary mask for current network

        if nnz(net_mask) == 0
            continue;  % Skip empty network masks
        end

        % Step 2a: Compute Dice score of each IC with current network
        dice_scores = zeros(num_ICs, 1);
        for i = 1:num_ICs
            ic_mask = binary_ICs(:,:,:,i);
            if nnz(ic_mask) == 0, continue; end

            intersection = nnz(ic_mask & net_mask);
            dice_scores(i) = 2 * intersection / (nnz(ic_mask) + nnz(net_mask));
        end

        % Step 2b: Sort ICs by Dice score (descending)
        [~, sorted_idx] = sort(dice_scores, 'descend');

        % Step 3: Greedily add ICs that improve Dice score
        selected_ICs = [];                        % ICs chosen for this network
        best_union_mask = false(size(net_mask)); % Union of selected IC masks
        best_dice = 0;                            % Best Dice score so far
        failures = 0;                             % Counter for consecutive non-improvements

        for k = 1:num_ICs
            this_ic = sorted_idx(k);
            ic_mask = binary_ICs(:,:,:,this_ic);

            % Combine current IC with existing union
            candidate_mask = best_union_mask | ic_mask;

            % Compute Dice score between candidate mask and network
            intersection = nnz(candidate_mask & net_mask);
            union_dice = 2 * intersection / (nnz(candidate_mask) + nnz(net_mask));

            if union_dice > best_dice
                % Accept this IC: it improves Dice
                best_dice = union_dice;
                best_union_mask = candidate_mask;
                selected_ICs(end+1) = this_ic;
                failures = 0;  % Reset failure counter
            else
                % IC did not improve overlap → count as failure
                failures = failures + 1;
            end

            % Stop if 3 ICs in a row don't help
            if failures >= 3
                break;
            end
        end

        % Step 4: Save results for this network
        final_ICs = [final_ICs; selected_ICs(:)];
        final_network_indices = [final_network_indices; repmat(net, numel(selected_ICs), 1)];
        final_network_labels = [final_network_labels; repmat(network_names(net), numel(selected_ICs), 1)];
    end

    % Step 5: Build final output table
    IC_assignment_table = table(final_ICs, final_network_indices, final_network_labels, ...
        'VariableNames', {'IC_Index', 'Assigned_Network', 'Network_Name'});

    % Optional: Display results in command window
    disp('Greedy Dice-based IC selection results:');
    disp(IC_assignment_table);
end
