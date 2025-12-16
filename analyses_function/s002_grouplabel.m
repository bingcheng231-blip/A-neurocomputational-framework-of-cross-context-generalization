function [group_label, rotation_direction, sub_traget_reward, exp_index] = s002_grouplabel(Re_analysis_data_path, exp_label)
%S002_GROUPLABEL 此处显示有关此函数的摘要
%   此处显示详细说明
    % VTC-dependent: group_label = 1;
    % random_distribution: group_label = 2;
    % VTC_depended_clockwise: rotation_direction = -1;
    % VTC_depended_anticlockwise: rotation_direction = 1;

    sub_reward = load([Re_analysis_data_path, 'd003_trials_reward.mat']).sub_reward;
    sub_info = load([Re_analysis_data_path, 'd004_Sub_learning_trials_info.mat']).sub_learning_trials_info;
    if strcmp(exp_label, 'all')
        sub_list = sub_info(:, 1);
        exp_index = true(36, 1);
    elseif strcmp(exp_label, 'MEG') || strcmp(exp_label, 'fMRI')
        exp_index = false(36, 1);
        for sub_i = 1:36
            exp_index(sub_i, 1) = contains(sub_info{sub_i, 1}, exp_label);
        end
        sub_list = sub_info(exp_index, 1);
    else
        error('exp_label: all, MEG, fMRI ')
    end


    cc_angle = (15:30:360)';
    cc_coord = [cosd(cc_angle), sind(cc_angle)];
    group_label = [];
    rotation_direction = zeros(size(sub_list, 1), 1);
    sub_traget_reward = cell(size(sub_list, 1), 1);
    for sub_i = 1:size(sub_list, 1)
        sub_name_list = strsplit(sub_list{sub_i, 1}, '_');
        sub_NO = str2num(sub_name_list{3});
        group_label(sub_i, 1) = mod(sub_NO, 2);
        trial_reward = sub_reward{sub_NO};
        sub_traget_reward{sub_i, 1} = trial_reward;
        if group_label(sub_i, 1) == 1
            if corr(trial_reward(:, 1), cc_coord(:, 1)) > 0 && corr(trial_reward(:, 2), cc_coord(:, 2))>0
                rotation_direction(sub_i, 1) = 1;
            elseif corr(trial_reward(:, 1), cc_coord(:, 1)) > 0 && corr(trial_reward(:, 2), cc_coord(:, 2))<0
                rotation_direction(sub_i, 1) = -1;
            elseif corr(trial_reward(:, 1), cc_coord(:, 1)) < 0 && corr(trial_reward(:, 2), cc_coord(:, 2))>0
                rotation_direction(sub_i, 1) = -1;
            elseif corr(trial_reward(:, 1), cc_coord(:, 1)) < 0 && corr(trial_reward(:, 2), cc_coord(:, 2))<0
                rotation_direction(sub_i, 1) = 1;
            end
        end
    end
    group_label(group_label==0)=2;

end

