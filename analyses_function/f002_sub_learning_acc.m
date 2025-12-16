function [outputArg1] = f002_sub_learning_acc(Re_analysis_data_path, exp_label)
%F002_SUB_LEARNING_ACC 此处显示有关此函数的摘要
%   此处显示详细说明

    sub_learning_trials_info = load([Re_analysis_data_path, 'd004_Sub_learning_trials_info.mat']).sub_learning_trials_info;
    [group_label, ~, ~, exp_index] = s002_grouplabel(Re_analysis_data_path, exp_label);
    sub_learning_trials = sub_learning_trials_info(exp_index, :);


    learning_session = zeros(15, size(sub_learning_trials, 1));
    learning_session_miss_trial = zeros(15, size(sub_learning_trials, 1));
    for sub_i = 1:size(sub_learning_trials, 1)
    
        run_acc = zeros(15, 1);
        miss_trial = zeros(15, 1);
        for run_i = 1:15
            stimuli_sequence = sub_learning_trials{sub_i, 2}{run_i};
            run_answer = sign(stimuli_sequence(:, 4));
            sub_report = sum(stimuli_sequence(:, 5:6), 2);
            miss_label = stimuli_sequence(:, 7)==1.5;
            out_rep_label = stimuli_sequence(:, 5)==0&stimuli_sequence(:, 6)==0;
            sub_report(miss_label) = 0;
            sub_report(out_rep_label) = 0;
    
            run_acc(run_i, 1) = sum(run_answer==sub_report)./72;
            miss_trial(run_i, 1) = sum(miss_label+out_rep_label);
        end
        learning_session(:, sub_i) = run_acc;
        learning_session_miss_trial(:, sub_i) = miss_trial;
    end
    
    % The original figure was generated using GraphPad.
    % VTC_depended_group
    g01_lr_acc = learning_session(:, group_label==1);
    % Random_distribution_group
    g02_lr_acc = learning_session(:, group_label==2);

    disp('Figure : Average acc of each run during the learning experiment')
    figure
    group_color = {[255, 128, 0]./256, [7, 126, 151]./256};
    y = g01_lr_acc'; x = 1:15;
    SEM=std(y)./sqrt(size(y, 1));
    s = shadedErrorBar(x,mean(y, 1),SEM, 'lineProps', '-');
    s.patch.FaceColor = [0.5, 0.5, 0.5];
    set(s.edge,'LineWidth',0.5, 'Color', [0.7, 0.7, 0.7])
    hold on
    plot(x, mean(y, 1), 'LineWidth', 2, 'color', group_color{1})
    hold on
    y = g02_lr_acc'; x = 1:15;
    SEM=std(y)./sqrt(size(y, 1));
    s = shadedErrorBar(x,mean(y, 1),SEM, 'lineProps', '-');
    s.patch.FaceColor = [0.5, 0.5, 0.5];
    set(s.edge,'LineWidth',0.5, 'Color', [0.7, 0.7, 0.7])
    hold on
    plot(x, mean(y, 1), 'LineWidth', 2, 'color', group_color{2})
    box off
    set(gca, 'FontSize', 15, 'FontWeight', 'bold')
    set(gca,'LineWidth', 2)

    outputArg1 = 'Analysis completed without error';
 end

