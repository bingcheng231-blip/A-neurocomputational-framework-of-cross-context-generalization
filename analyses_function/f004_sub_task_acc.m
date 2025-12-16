function [outputArg1] = f004_sub_task_acc(Re_analysis_data_path, exp_label)
%F004_SUB_TASK_ACC 此处显示有关此函数的摘要
%   此处显示详细说明

    sub_task_trials_info = load([Re_analysis_data_path, 'd005_Sub_task_trials_info.mat']).sub_task_trials_info;
    [group_label, ~, ~, exp_index] = s002_grouplabel(Re_analysis_data_path, exp_label);
    sub_task_trials = sub_task_trials_info(exp_index, :);

    task_session = zeros(10, size(sub_task_trials, 1));
    task_session_miss_trial = zeros(10, size(sub_task_trials, 1));
    for sub_i = 1:size(sub_task_trials, 1)
    
        run_acc = zeros(10, 1);
        miss_trial = zeros(10, 1);
        for run_i = 1:10
            stimuli_sequence = sub_task_trials{sub_i, 2}{run_i};
            run_answer = sign(stimuli_sequence(:, 4));
            sub_report = sum(stimuli_sequence(:, 5:6), 2);
            miss_label = stimuli_sequence(:, 7)==1.5;
            out_rep_label = stimuli_sequence(:, 5)==0&stimuli_sequence(:, 6)==0;
            sub_report(miss_label) = 0;
            sub_report(out_rep_label) = 0;
    
            run_acc(run_i, 1) = sum(run_answer==sub_report)./72;
            miss_trial(run_i, 1) = sum(miss_label+out_rep_label);
        end
        task_session(:, sub_i) = run_acc;
        task_session_miss_trial(:, sub_i) = miss_trial;
    end
    
    % The original figure was generated using GraphPad.
    % VTC_dependent_group
    group_acc{1} = mean(task_session(:, group_label==1));
    % Random_distribution_group
    group_acc{2} = mean(task_session(:, group_label==2));

    [p,h,stats] = ranksum(group_acc{1}, group_acc{2}, 'tail', 'right');
    disp(['task acc, ranksum, VTC dependent > Random distribution: pvalue=', num2str(p)])

    % Average acc of the task experiment
    disp(' ')
    disp('Figure : Average acc of the task experiment')
    group_name = {'VTC-dependent', 'random distribution'};
    group_color = {[255, 128, 0]./256, [7, 126, 151]./256};
    figure;
    for gl = 1:2
        group_data = group_acc{gl};
        group_means = mean(group_data);
        group_sems  = std(group_data)/sqrt(size(group_data, 1));
        jitter = 0.2;
        b = bar(gl-0.5, group_means, 0.6, 'FaceColor','none', 'EdgeColor','none');
        hold on
        xi = gl-0.5 + (rand(size(group_data)) - 0.5) * 2 * jitter;
        scatter(xi, group_data, 50, ones(1, 3).*0.2, 'filled', ...
            'MarkerFaceAlpha',0.7, 'MarkerEdgeAlpha',0.2);
        b.EdgeColor = group_color{gl};       
        b.LineWidth = 5;
        errorbar(gl-0.5, group_means, group_sems, 'color', group_color{gl}, 'LineWidth',3, 'CapSize',20, 'LineStyle','none');
    end
    box off
    set(gca, 'FontSize', 15, 'FontWeight', 'bold')
    set(gca,'LineWidth', 2)
    ylim([0 1])
    xticks(0.5:3);
    xticklabels(group_name);
    disp(' ')
    outputArg1 = 'Analysis completed without error';

end

