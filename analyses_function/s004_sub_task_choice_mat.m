function [sub_task_choice_mat] = s004_sub_task_choice_mat(Re_analysis_data_path,exp_label)
%S004_SUB_CHOICE_MAT 此处显示有关此函数的摘要
%   此处显示详细说明
    sub_task_trials_info = load([Re_analysis_data_path, 'd005_Sub_task_trials_info.mat']).sub_task_trials_info;
    [~, ~, ~, exp_index] = s002_grouplabel(Re_analysis_data_path, exp_label);
    sub_task_trials = sub_task_trials_info(exp_index, :);

    sub_task_choice_mat = cell(size(sub_task_trials, 1), 1);
    for sub_i = 1:size(sub_task_trials, 1)
        % [theta, ctxA_reward, ctxB_reward, bad_hit_A, bad_hit_B, good_hit_A, good_hit_B, resp_t_A, resp_t_B]   
        sub_choice_mat = zeros(12, 9);
        for run_i = 1:10
            stimuli_sequence = sub_task_trials{sub_i, 2}{run_i};
            for center_i = 1:12
                % center label
                center_label = stimuli_sequence(:, 3) == center_i;
                center_choice = [stimuli_sequence(center_label, [2, 4:7])];

                % center coordinate in reward space
                center_coordinate_a = unique(center_choice(center_choice(:, 1) == 1, 2));
                center_coordinate_b = unique(center_choice(center_choice(:, 1) == 2, 2));
                [theta,rho] = cart2pol (center_coordinate_a,center_coordinate_b);

                good_hit_A = sum(center_choice(center_choice(:, 1) == 1, 3));
                good_hit_B = sum(center_choice(center_choice(:, 1) == 2, 3));

                bad_hit_A = sum(center_choice(center_choice(:, 1) == 1, 4));
                bad_hit_B = sum(center_choice(center_choice(:, 1) == 2, 4));

                sub_choice_mat(center_i, 1) = theta;
                sub_choice_mat(center_i, 2) = center_coordinate_a;
                sub_choice_mat(center_i, 3) = center_coordinate_b;

                sub_choice_mat(center_i, 4) = sub_choice_mat(center_i, 4) + bad_hit_A;
                sub_choice_mat(center_i, 5) = sub_choice_mat(center_i, 5) + bad_hit_B;

                sub_choice_mat(center_i, 6) = sub_choice_mat(center_i, 6) + good_hit_A;
                sub_choice_mat(center_i, 7) = sub_choice_mat(center_i, 7) + good_hit_B;

                run_miss_label = center_choice(:, 5)>=1.5;

                resp_t_A = mean(center_choice(center_choice(:, 1) == 1 & ~run_miss_label, 5));
                resp_t_B = mean(center_choice(center_choice(:, 1) == 2 & ~run_miss_label, 5));

                sub_choice_mat(center_i, 8) = mean([sub_choice_mat(center_i, 8), resp_t_A]);
                sub_choice_mat(center_i, 9) = mean([sub_choice_mat(center_i, 9), resp_t_B]);
            end
        end
        sub_choice_mat(:, 4:7) = sub_choice_mat(:, 4:7)./(3*run_i);
        sub_task_choice_mat{sub_i, 1} = sub_choice_mat;
    end
end

