function IC_assignment_table = assign_ICs_to_networks_clust(networks, IC_masks, network_names)
% assign_ICs_to_networks - Assigns ICs to networks using k-means on Dice scores
%
% Inputs:
%   networks       - 4D binary matrix (X x Y x Z x N), N = number of networks
%   IC_masks       - 4D binary matrix (X x Y x Z x M), M = number of ICs
%   network_names  - Cell array of strings with names for each network (length N)
%
% Output:
%   IC_assignment_table - Table with selected ICs, assigned network index, and name

    num_networks = size(networks, 4);
    num_ICs = size(IC_masks, 4);

    % Validate network_names
    if nargin < 3 || isempty(network_names)
        network_names = strcat("Network_", string(1:num_networks));
    elseif length(network_names) ~= num_networks
        error('Length of network_names must match number of networks (4th dim of `networks`).');
    end

    dice_matrix = zeros(num_ICs, num_networks);

    % Step 1: Compute Dice coefficients
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

    % Step 2: Assign each IC to the best network (highest Dice)
    [best_dice, best_net_idx] = max(dice_matrix, [], 2);
    IC_assignment_mask = false(num_ICs, num_networks);
    for ic = 1:num_ICs
        if best_dice(ic) > 0
            IC_assignment_mask(ic, best_net_idx(ic)) = true;
        end
    end

    % Step 3: For each network, use k-means clustering and select strongest group
    final_ICs = [];
    final_network_indices = [];
    final_network_labels = [];

    for net = 1:num_networks
        assigned_ICs = find(IC_assignment_mask(:, net));
        if isempty(assigned_ICs), continue; end

        dice_scores = dice_matrix(assigned_ICs, net);

        if length(dice_scores) < 3
            [~, idx] = max(dice_scores);
            selected = assigned_ICs(idx);
        else
            sorted_dice = sort(dice_scores);
            init_centroids = [sorted_dice(1); median(dice_scores); sorted_dice(end)];

            [cluster_idx, centroids] = kmeans(dice_scores, 3, 'Start', init_centroids);
            [~, strongest_cluster] = max(centroids);
            selected = assigned_ICs(cluster_idx == strongest_cluster);
        end

        final_ICs = [final_ICs; selected(:)];
        final_network_indices = [final_network_indices; repmat(net, length(selected), 1)];
        final_network_labels = [final_network_labels; repmat(network_names(net), length(selected), 1)];
    end

    % Step 4: Output table
    IC_assignment_table = table(final_ICs, final_network_indices, final_network_labels, ...
        'VariableNames', {'IC_Index', 'Assigned_Network', 'Network_Name'});

    disp('Final selected ICs (using k-means) and their assigned networks:');
    disp(IC_assignment_table);
end
