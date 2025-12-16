function [outputArg1,outputArg2] = f007_vmPFC_feature_rotation_and_compression(Re_analysis_data_path, exp_label)
%F007_VMPFC_FEATURE_ROTATION_AND_COMPRESSION 此处显示有关此函数的摘要
%   此处显示详细说明
    [group_label, ~, ~, exp_index] = s002_grouplabel(Re_analysis_data_path, exp_label);
    Sub_fMRI_ROI_RDM = load([Re_analysis_data_path, 'd006_Sub_fMRI_ROIs_RDM.mat']);
    % Sub_fMRI_ROI_RDM.roi_rdm: sub_num * ROI_num * condition(day01, day03ctxA, day03ctxB, day03all)
    % subject list => Sub_fMRI_ROI_RDM.sub_list
    % ROI list => Sub_fMRI_ROI_RDM.roi_name
    % for each rdm, row 1~12 => category 1~12; row 1~24 => [ctxA category 1~12; ctxB category 1~12]
    % Due to damaged data labels, the fMRI data of 'Sub_02_02_MEG_fMRI'  on the day 01 is missing
    STRUCT_flag = true;
    all_behavior_model = load([Re_analysis_data_path, 'd007_behavior_model.mat']).all_behavior_model;
    [sub_structure_model_objectEV] = s005_expected_value_of_object(Re_analysis_data_path, exp_label, STRUCT_flag, all_behavior_model);

    roi_cc_rdm = Sub_fMRI_ROI_RDM.roi_rdm;
    roi_num = 2; % vmPFC
    sub_feature_list = {};
    for sub_i = 1:size(sub_structure_model_objectEV, 1)
        structure_model_EV(1, :) = sub_structure_model_objectEV{sub_i, 1}(end, :);
        sub_feature_list{sub_i, 1} = [structure_model_EV(1, 1:12)', structure_model_EV(1, 13:24)'];
    end
    [roi_features_dist, roi_rotation_rr] = s007_rdm_feature_rotation_and_compression(roi_cc_rdm, roi_num, sub_feature_list);


      
    disp(' ')
    disp('figure : Geometric relationships of category’s EVs across different contexts in the vmPFC on Day 3')
    figure
    theta = deg2rad(1:361);
    group_color = {[255, 128, 0]./256, [7, 126, 151]./256};
    ShadowAlpha = 0.2;
    background_MaxCorr = 0.8;
    background_LineWidth = 1.2;
    background_Linecolor = 0.7.*ones(1, 3);
    for ri = 0.2:0.2:background_MaxCorr
        [x,y] = pol2cart(theta, ri);
        plot(x, y, 'LineWidth', background_LineWidth, 'color',  background_Linecolor)
        hold on
    end
    for ti = 0:30:330
        [x,y] = pol2cart(deg2rad(ti), (0:0.01:background_MaxCorr)');
        plot(x, y, 'LineWidth', background_LineWidth, 'color',  background_Linecolor)
        hold on
    end
    for gl = 1:2
        rotation_rr = squeeze(roi_rotation_rr(group_label==gl, :));
        SEM=std(rotation_rr)./sqrt(sum(group_label==gl));
        [x1_SEM,y1_SEM] = pol2cart(theta, mean(rotation_rr)-SEM);
        [x2_SEM,y2_SEM] = pol2cart(theta, mean(rotation_rr)+SEM);
        fill([x1_SEM, fliplr(x2_SEM)], [y1_SEM, fliplr(y2_SEM)], 'k', 'FaceAlpha', ShadowAlpha, 'EdgeColor', 'none');
        hold on
        [x,y] = pol2cart(theta, mean(rotation_rr));
        plot(x, y, 'LineWidth', 1.7, 'color', group_color{gl})
        hold on
    end
    axis equal
    axis off

    % MDS
    disp(' ')
    disp('Visualization of vmPFCs representation of EVs in a 3D space using MDS')
    rearrange_rdm = [];
    for sub_i = 1:27
        structure_model_EV(1, :) = sub_structure_model_objectEV{sub_i, 1}(end, :);
        value_sequence = [sortrows([(1:12)', structure_model_EV(1, 1:12)',], 2); sortrows([(13:24)', structure_model_EV(1, 13:24)'], 2)];
        roi_rdm = squareform(roi_cc_rdm{sub_i, roi_num, 4});
        rearrange_rdm(:, :, sub_i) = roi_rdm(value_sequence(:, 1), value_sequence(:, 1));
    end
    value_colormap = jet(360);
    value_color = value_colormap(1:30:360, :);
    group_name = {'VTC_depended', 'Random_distribution'};
    for gl = 1:2
        disp(['figure : ', group_name{gl}])
        group_average_rdm = squeeze(mean(rearrange_rdm(:, :, group_label==gl), 3));
        [Y,e] = cmdscale(group_average_rdm, 3);
        figure
        scatter3(Y(1:12, 1), Y(1:12, 2), Y(1:12, 3), 100, value_color, 'filled', "o"	)
        hold on
        scatter3(Y(13:24, 1), Y(13:24, 2), Y(13:24, 3), 100, value_color, 'filled', "^")
        box off
        set(gca, 'FontSize', 15, 'FontWeight', 'bold')
        set(gca,'LineWidth',2)
    end

    sub_roi_ra = [];
    sa_theta = deg2rad(-170:10:180)';
    for sub_i = 1:27
        ra = roi_rotation_rr(sub_i, 10:10:360)';
        sub_ra = [ra(19:36); ra(1:18)];
        [x_list] = s009_Gaussian_perfer_direction(sub_ra,  [] , sa_theta, 32);
        sub_roi_ra(sub_i, 1) = rad2deg(x_list(1, 3));
    end
    disp(' ')
    disp('Group level cross-context EVs rotation in the vmPFC:')
    for gl = 1:2
        disp(group_name{gl})
        group_ra = sub_roi_ra(group_label==gl, 1);
        mean_rotation = mean(group_ra);
        SEM_rotation=std(group_ra)./sqrt(size(group_ra, 1));
        disp(['mean_rotation: ', num2str(mean_rotation)]);
        disp(['sem: ', num2str(SEM_rotation)])
        disp(' ')
    end


    roi_num = 2; % vmPFC
    cc_angle = (15:30:360)';
    cc_coord = [cosd(cc_angle), sind(cc_angle)];
    sub_feature_list = {};
    for sub_i = 1:size(sub_structure_model_objectEV, 1)
        sub_feature_list{sub_i, 1} = cc_coord;
    end
    [roi_features_dist, roi_rotation_rr] = s007_rdm_feature_rotation_and_compression(roi_cc_rdm, roi_num, sub_feature_list);

    disp(' ')
    disp('figure : Geometric relationships of VTC features across different contexts in the vmPFC on Day 3')
    figure
    theta = deg2rad(1:361);
    group_color = {[255, 128, 0]./256, [7, 126, 151]./256};
    ShadowAlpha = 0.2;
    background_MaxCorr = 0.8;
    background_LineWidth = 1.2;
    background_Linecolor = 0.7.*ones(1, 3);
    for ri = 0.2:0.2:background_MaxCorr
        [x,y] = pol2cart(theta, ri);
        plot(x, y, 'LineWidth', background_LineWidth, 'color',  background_Linecolor)
        hold on
    end
    for ti = 0:30:330
        [x,y] = pol2cart(deg2rad(ti), (0:0.01:background_MaxCorr)');
        plot(x, y, 'LineWidth', background_LineWidth, 'color',  background_Linecolor)
        hold on
    end
    for gl = 1:2
        rotation_rr = squeeze(roi_rotation_rr(group_label==gl, :));
        SEM=std(rotation_rr)./sqrt(sum(group_label==gl));
        [x1_SEM,y1_SEM] = pol2cart(theta, mean(rotation_rr)-SEM);
        [x2_SEM,y2_SEM] = pol2cart(theta, mean(rotation_rr)+SEM);
        fill([x1_SEM, fliplr(x2_SEM)], [y1_SEM, fliplr(y2_SEM)], 'k', 'FaceAlpha', ShadowAlpha, 'EdgeColor', 'none');
        hold on
        [x,y] = pol2cart(theta, mean(rotation_rr));
        plot(x, y, 'LineWidth', 1.7, 'color', group_color{gl})
        hold on
    end
    axis equal
    axis off

    sub_roi_ra = [];
    sa_theta = deg2rad(-170:10:180)';
    for sub_i = 1:27
        ra = roi_rotation_rr(sub_i, 10:10:360)';
        sub_ra = [ra(19:36); ra(1:18)];
        [x_list] = s009_Gaussian_perfer_direction(sub_ra,  [] , sa_theta, 32);
        sub_roi_ra(sub_i, 1) = rad2deg(x_list(1, 3));
    end
    disp(' ')
    disp('Group level cross-context VTC features rotation in the vmPFC:')
    for gl = 1:2
        disp(group_name{gl})
        group_ra = sub_roi_ra(group_label==gl, 1);
        mean_rotation = mean(group_ra);
        SEM_rotation=std(group_ra)./sqrt(size(group_ra, 1));
        disp(['mean_rotation: ', num2str(mean_rotation)]);
        disp(['sem: ', num2str(SEM_rotation)])
        disp(' ')
    end

end

