function [sub_objectEV] = s005_expected_value_of_object(Re_analysis_data_path, exp_label, STRUCT_flag, all_behavior_model)
%S005_EXPECTED_VALUE_OF_OBJECT 此处显示有关此函数的摘要
%   此处显示详细说明
    sub_learning_trials_info = load([Re_analysis_data_path, 'd004_Sub_learning_trials_info.mat']).sub_learning_trials_info;
    [~, ~, ~, exp_index] = s002_grouplabel(Re_analysis_data_path, exp_label);
    sub_learning_trials = sub_learning_trials_info(exp_index, :);
    all_behavior_model = all_behavior_model(exp_index, :);

    sub_objectEV = cell(size(sub_learning_trials, 1), 1);
    decayed_lr_flag = true;
    prior_structure = [];
    for sub_i = 1:size(sub_learning_trials, 1)
        sub_choices = cell2mat(sub_learning_trials{sub_i, 2}');
        stimuli_id = (sub_choices(:, 2)-1)*12+sub_choices(:, 3);
        outcomes = sub_choices(:, 4);
        choices = sub_choices(:, 5);
        if STRUCT_flag
            model_theta = all_behavior_model{sub_i, 2};
        else
            model_theta = all_behavior_model{sub_i, 1};
        end
        [negLogLik,VV,PP, step_lr] = s003_ctx_RL_LogLikelihood(model_theta(1, 1:end-1),stimuli_id,outcomes,choices,STRUCT_flag, decayed_lr_flag, prior_structure);
        sub_objectEV{sub_i, 1} = VV;
        sub_objectEV{sub_i, 2} = PP;
    end
end

