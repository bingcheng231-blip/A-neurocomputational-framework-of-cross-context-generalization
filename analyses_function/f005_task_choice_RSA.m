function [outputArg1] = f005_task_choice_RSA(Re_analysis_data_path,exp_label)
%F005_TASK_CHOICE_RSA 此处显示有关此函数的摘要
%   此处显示详细说明
    [sub_task_choice_mat] = s004_sub_task_choice_mat(Re_analysis_data_path,exp_label);
    [group_label, ~, ~, exp_index] = s002_grouplabel(Re_analysis_data_path, exp_label);
    sub_task_choice_mat = sub_task_choice_mat(exp_index, :);
    % day03 subject choice RSA fitting
    sub_beta = nan(size(sub_task_choice_mat, 1), 2);
    for sub_i = 1:size(sub_task_choice_mat, 1)
        sub_choice_mat = sub_task_choice_mat{sub_i, 1};

        sub_choice_mat = sortrows(sub_choice_mat, 1);
        B_model_A = double(sub_choice_mat(:, 2)>0); % Good in context A
        B_model_B = double(sub_choice_mat(:, 3)>0); % Good in context B
        B_model = [B_model_A; B_model_B];
        boundary_RDM = zscore(pdist(B_model))';

        L_model = zeros(12, 2);
        L_model(sub_choice_mat(:, 2)>-sub_choice_mat(:, 3), 1) = 1;
        L_model(sub_choice_mat(:, 2)<-sub_choice_mat(:, 3), 2) = 1;
        L_model(sub_choice_mat(:, 2)==-sub_choice_mat(:, 3), 1:2) = 1;
        L_model_AB = [L_model;L_model];
        linear_model = zscore(pdist(L_model_AB))';

        % day 03
        hit_conext_A = sub_choice_mat(:, 6);
        hit_conext_B = sub_choice_mat(:, 7);
        choice_matrix = [hit_conext_A; hit_conext_B];
        day_03_choice_RDM = zscore(pdist(choice_matrix))';
        b = regress(day_03_choice_RDM, [boundary_RDM, linear_model]);
        sub_beta(sub_i, :) = b';
    end


    disp(['figure : task choice RSA'])
    disp(' ')
    figure;
    group_color = {[255, 128, 0]./256, [7, 126, 151]./256};
    group_beta = sub_beta(group_label==1, :);
    figure;
    gl = 1;
    for mi = 1:size(group_beta, 2)
        gm_beta = group_beta(:, mi);
        gm_means = mean(gm_beta);
        gm_sems  = std(gm_beta)/sqrt(size(gm_beta, 1));
        jitter = 0.2;
        b = bar(mi-0.5, gm_means, 0.6, 'FaceColor','none', 'EdgeColor','none');
        hold on
        xi = mi-0.5 + (rand(size(gm_beta)) - 0.5) * 2 * jitter;
        scatter(xi, gm_beta, 50, ones(1, 3).*0.2, 'filled', ...
            'MarkerFaceAlpha',0.7, 'MarkerEdgeAlpha',0.2);
        b.EdgeColor = group_color{gl};
        b.LineWidth = 5;
        errorbar(mi-0.5, gm_means, gm_sems, 'color', group_color{gl}, 'LineWidth',2, 'CapSize',12, 'LineStyle','none');
    end
    gl = 2;
    group_beta = sub_beta(group_label==gl, :);
    for mi = 1:size(group_beta, 2)
        gm_beta = group_beta(:, mi);
        gm_means = mean(gm_beta);
        gm_sems  = std(gm_beta)/sqrt(size(gm_beta, 1));
        jitter = 0.2;
        b = bar(mi-0.5+2, gm_means, 0.6, 'FaceColor','none', 'EdgeColor','none');
        hold on
        xi = mi-0.5+2 + (rand(size(gm_beta)) - 0.5) * 2 * jitter;
        scatter(xi, gm_beta, 50, ones(1, 3).*0.2, 'filled', ...
            'MarkerFaceAlpha',0.7, 'MarkerEdgeAlpha',0.2);
        b.EdgeColor = group_color{gl};
        b.LineWidth = 5;
        errorbar(mi-0.5+2, gm_means, gm_sems, 'color', group_color{gl}, 'LineWidth',2, 'CapSize',12, 'LineStyle','none');
    end
    box off
    set(gca, 'FontSize', 15, 'FontWeight', 'bold')
    set(gca,'LineWidth', 2)
    ylim([-0.5 1])
    xticks(0.5:4);
    xticklabels({'Boundary', 'Linear', 'Boundary', 'Linear'});

    outputArg1 = 'Analysis completed without error';
end

