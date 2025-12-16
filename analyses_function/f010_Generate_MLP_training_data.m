function [outputArg1] = f010_Generate_MLP_training_data(Re_analysis_data_path,exp_label)
%F010_GENERATE_MLP_TRAINING_DATA 此处显示有关此函数的摘要
%   此处显示详细说明
    [~, ~, sub_traget_reward, ~] = s002_grouplabel(Re_analysis_data_path, exp_label);
    % Sub_fMRI_ROI_RDM.roi_rdm: sub_num * ROI_num * condition(day01, day03ctxA, day03ctxB, day03all)
    % subject list => Sub_fMRI_ROI_RDM.sub_list
    % ROI list => Sub_fMRI_ROI_RDM.roi_name
    % for each rdm, row 1~12 => category 1~12; row 1~24 => [ctxA category 1~12; ctxB category 1~12]
    % Due to damaged data labels, the fMRI data of 'Sub_02_02_MEG_fMRI'  on the day 01 is missing
    STRUCT_flag = true;
    all_behavior_model = load([Re_analysis_data_path, 'd007_behavior_model.mat']).all_behavior_model;
    [sub_structure_model_objectEV] = s005_expected_value_of_object(Re_analysis_data_path, exp_label, STRUCT_flag, all_behavior_model);


    cc_angle = (15:30:360)';
    cc_coord = [cosd(cc_angle), sind(cc_angle)];
    sub_data = cell(size(sub_structure_model_objectEV, 1), 2);
    for sub_i = 1:size(sub_structure_model_objectEV, 1)
        trial_reward = [sub_traget_reward{sub_i}(:, 1); sub_traget_reward{sub_i}(:, 2)];
        sub_EV = sub_structure_model_objectEV{sub_i, 1}(end, :)';
        test_trial_xy = repmat([cosd(cc_angle), sind(cc_angle)], 2, 1);
        test_trial = [test_trial_xy, [[ones(12, 1), zeros(12, 1)];[zeros(12, 1), ones(12, 1)]]];
    
        ctxA_cc = [cc_coord, ones(12, 1), zeros(12, 1)];
        ctxB_cc = [cc_coord, zeros(12, 1), ones(12, 1)];
        ctx_cc = [ctxA_cc; ctxB_cc];
        ctx_cc = [ctx_cc, trial_reward, sub_EV];
    
        vtc_cc = repmat(ctx_cc, 4000, 1);
        vtc_cc(:, 1:2) = vtc_cc(:, 1:2)+0.4.*rand(size(vtc_cc, 1), 2)-0.2;
        vtc_cc = vtc_cc(randperm(size(vtc_cc, 1)), :);
    
        sub_data{sub_i, 1} = vtc_cc(:, 1:4); % original cc
    
        sub_data{sub_i, 3} = vtc_cc(:, 5); % target rewards
        sub_data{sub_i, 4} = vtc_cc(:, 6); % EVs from the structure model
    end
    % save the sub_data for the MLP training and the test trial for the MLP testing

    disp(' ')
    outputArg1 = 'Analysis completed without error';
end

