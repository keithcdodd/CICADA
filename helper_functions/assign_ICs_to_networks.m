function IC_assignment_table = assign_ICs_to_networks(networks, IC_masks)
% assign_ICs_to_networks - Assigns ICs to networks using greedy Dice improvement strategy
% after restricting ICs to only the network they best match individually.
%
% Inputs:
%   networks   - 4D binary matrix (X x Y x Z x N), N = number of networks
%   IC_masks   - 4D binary matrix (X x Y x Z x M), M = number of ICs
%
% Output:
%   IC_assignment_table - Table with final selected ICs and their assigned network

    num_networks = size(networks, 4);
    num_ICs = size(IC_masks, 4);
    dice_matrix = zeros(num_ICs, num_networks);

    % Step 1: Compute individual Dice coefficients
    for ic = 1:num_ICs
        IC_mask = logical(IC_masks(:,:,:,ic));
        IC_voxels = sum(IC_mask(:));
        if IC_voxels == 0, continue; end

        for net = 1:num_networks
            net_mask = logical(networks(:,:,:,net));
            net_voxels = sum(net_mask(:));
            if net_voxels == 0, continue; end

            intersection = sum(IC_mask & net_mask, 'all');
            dice_matrix(ic, net) = 2 * intersection / (IC_voxels + net_voxels);
        end
    end

    % Step 2: Assign each IC to the network with highest individual Dice
    [best_dice, best_net_idx] = max(dice_matrix, [], 2);
    IC_assignment_mask = false(num_ICs, num_networks);
    for ic = 1:num_ICs
        if best_dice(ic) > 0
            IC_assignment_mask(ic, best_net_idx(ic)) = true;
        end
    end

    % Step 3: For each network, perform cumulative selection among assigned ICs
    final_ICs = [];
    final_networks = [];

    for net = 1:num_networks
        curr_network = logical(networks(:,:,:,net));
        assigned_ICs = find(IC_assignment_mask(:, net));

        if isempty(assigned_ICs)
            continue;
        end

        % Sort assigned ICs by individual Dice with this network
        [~, sort_idx] = sort(dice_matrix(assigned_ICs, net), 'descend');
        sorted_ICs = assigned_ICs(sort_idx);

        % Initialize cumulative mask and Dice
        cum_mask = false(size(curr_network));
        prev_dice = 0;
        selected_ICs = [];

        for i = 1:length(sorted_ICs)
            ic = sorted_ICs(i);
            new_mask = cum_mask | logical(IC_masks(:,:,:,ic));
            intersection = sum(new_mask & curr_network, 'all');
            union = sum(new_mask | curr_network, 'all');
            new_dice = intersection / union;

            if new_dice > prev_dice
                cum_mask = new_mask;
                prev_dice = new_dice;
                selected_ICs(end+1) = ic;
            else
                break; % Stop if Dice does not improve
            end
        end

        % Store selected ICs and their network
        final_ICs = [final_ICs; selected_ICs(:)];
        final_networks = [final_networks; repmat(net, length(selected_ICs), 1)];
    end

    % Step 4: Output final ICs and assigned networks as a table
    IC_assignment_table = table(final_ICs, final_networks, ...
        'VariableNames', {'IC_Index', 'Assigned_Network'});

    disp('Final selected ICs and their assigned networks:');
    disp(IC_assignment_table);
end
