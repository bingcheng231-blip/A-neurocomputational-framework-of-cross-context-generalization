function [outputArg1,outputArg2] = f009_VTC_feature_rotation_and_compression(Re_analysis_data_path,exp_label,figure_generation_data_path, re_analysis)

%F009_VTC_FEATURE_ROTATION_AND_COMPRESSION 此处显示有关此函数的摘要
%   此处显示详细说明
    [group_label, rotation_direction, ~, exp_index] = s002_grouplabel(Re_analysis_data_path, exp_label);
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
    cc_angle = (15:30:360)';
    cc_coord = [cosd(cc_angle), sind(cc_angle)];
    VTC_cc_rdm(:, 1) = zscore(pdist(cc_coord(:, 1)))';
    VTC_cc_rdm(:, 2) = zscore(pdist(cc_coord(:, 2)))';
    
    
    ctx_cc_beta = nan(27, 3, 3);
    roi_num = 5; % VTC
    for sub_i = 1:27
        for ctx_i = 1:3
            roi_rdm = zscore(roi_cc_rdm{sub_i, roi_num, ctx_i}');
            if ~isempty(roi_rdm)
                ctx_cc_beta(sub_i, ctx_i, :) = regress(roi_rdm, [VTC_cc_rdm, ones(66, 1)]);
            end
        end
    end

    disp(['figure : Representations of VTC feature 1 and VTC feature 2 under different conditions in the VTC. '])
    disp(' ')
    context_color{1} = [0.4, 0.4, 0.4]; context_color{2} = [1, 0.15, 0.85]; context_color{3} = [0.15, 0.85, 1];
    group_name = {'VTC_depended_group', 'Random_distribution_group'};
    bar_name = {'day01 VTC f01', 'day01 VTC f02', 'day03 CtxA VTC f01', 'day03 CtxA VTC f02', ...
        'day03 CtxB VTC f01', 'day03 CtxB VTC f02'};
    for gl = 1:2
        disp([group_name{gl}, ':'])
        bn = 1;
        figure
        for ctx_i = 1:3
            for fi = 1:2
                gm_beta = ctx_cc_beta(group_label==gl, ctx_i, fi);
                gm_beta = gm_beta(~isnan(gm_beta));
                gm_means = mean(gm_beta);
                gm_sems  = std(gm_beta)/sqrt(size(gm_beta, 1));
                % jitter = 0.2;
                b = bar(bn-0.5, gm_means, 0.6, 'FaceColor','none', 'EdgeColor','none');
                hold on
                % xi = bn-0.5 + (rand(size(gm_beta)) - 0.5) * 2 * jitter;
                % scatter(xi, gm_beta, 50, ones(1, 3).*0.2, 'filled', ...
                %     'MarkerFaceAlpha',0.7, 'MarkerEdgeAlpha',0.2);
                b.EdgeColor = context_color{ctx_i};
                b.LineWidth = 5;
                errorbar(bn-0.5, gm_means, gm_sems, 'color', context_color{ctx_i}, 'LineWidth',2, 'CapSize',12, 'LineStyle','none');
                bn = bn+1;
            end
        end
        box off
        set(gca, 'FontSize', 15, 'FontWeight', 'bold')
        set(gca,'LineWidth', 2)
        if gl == 1
            ylim([0 1])
        else
            ylim([0 0.8])
        end
        xticks(1:2:14);
        xticklabels(bar_name);
    end
 
    for gl = 1:2
        disp([group_name{gl}, ':'])
        for ctx_i = 1:3
            if ctx_i == 3
                [h,p,ci,stats] = ttest(ctx_cc_beta(group_label==gl, ctx_i, 1), ctx_cc_beta(group_label==gl, ctx_i, 2), 'tail', 'left');
                disp(['ttest, ', bar_name{ctx_i*2-1}, '<', bar_name{ctx_i*2}, ': p-value=', num2str(p)])
            else
                [h,p,ci,stats] = ttest(ctx_cc_beta(group_label==gl, ctx_i, 1), ctx_cc_beta(group_label==gl, ctx_i, 2), 'tail', 'right');
                disp(['ttest, ', bar_name{ctx_i*2-1}, '>', bar_name{ctx_i*2}, ': p-value=', num2str(p)])
                tvalue(gl, ctx_i) = stats.tstat;
            end
        end
        disp(' ')
    end



    % VTC feature model comparison through RSA
    roi_num = 5;
    sub_roi_beta = [];
    for sub_i = 1:27
        cc_grid = [cc_coord; cc_coord];
        cc_orth = [[cc_coord(:, 1), zeros(12, 1)]; [zeros(12, 1), cc_coord(:, 2)]];
        if rotation_direction(sub_i, 1) == -1
            cc_para = [[cc_coord(:, 1), zeros(12, 1)]; [-cc_coord(:, 2), zeros(12, 1)]];
        else
            cc_para = [[cc_coord(:, 1), zeros(12, 1)]; [cc_coord(:, 2), zeros(12, 1)]];
        end
        cc_rdm(:, 1) = zscore(pdist(cc_grid))';
        cc_rdm(:, 2) = zscore(pdist(cc_orth))';
        cc_rdm(:, 3) = zscore(pdist(cc_para))';
    
        roi_rdm = zscore(roi_cc_rdm{sub_i, roi_num, 4}');
        [b,bint,r] = regress(roi_rdm, [cc_rdm, ones(276, 1)]);
        sub_roi_beta(:, sub_i) = b;
    end

    disp(['figure : fitting the VTC feature with different structures to the VTC'])
    disp(' ')
    figure;
    group_color = {[255, 128, 0]./256, [7, 126, 151]./256};
    VTC_subgroup_color = {[255, 216, 178]./256, [255, 164, 164]./256};
    three_group_color = [VTC_subgroup_color, group_color(2)];
    group_name = {'VTC_depended_group', 'Random_distribution_group'};
    VTC_subgroup_name = {'VTC_depended_clockwise', 'VTC_depended_anticlockwise'};
    three_group_name = [VTC_subgroup_name, group_name(2)];
    bn = 1;
    bar_name = {'Grid', 'Orthogonal', 'Parallel'};
    rotation_label = rotation_direction;
    rotation_label(rotation_label==0) = 3;
    rotation_label(rotation_label==1) = 2;
    rotation_label(rotation_label==-1) = 1;
    for bi = 1:3
        for gl = 1:3
            gm_beta = sub_roi_beta(bi, rotation_label==gl);
            gm_beta = gm_beta(~isnan(gm_beta));
            gm_means = mean(gm_beta);
            gm_sems  = std(gm_beta)/sqrt(size(gm_beta, 2));
            % jitter = 0.2;
            b = bar(bn-0.5, gm_means, 0.6, 'FaceColor','none', 'EdgeColor','none');
            hold on
            % xi = bn-0.5 + (rand(size(gm_beta)) - 0.5) * 2 * jitter;
            % scatter(xi, gm_beta, 50, ones(1, 3).*0.2, 'filled', ...
            %     'MarkerFaceAlpha',0.7, 'MarkerEdgeAlpha',0.2);
            b.EdgeColor = three_group_color{gl};
            b.LineWidth = 5;
            errorbar(bn-0.5, gm_means, gm_sems, 'color', three_group_color{gl}, 'LineWidth',2, 'CapSize',12, 'LineStyle','none');
            bn = bn+1;
        end
    end
    box off
    set(gca, 'FontSize', 15, 'FontWeight', 'bold')
    set(gca,'LineWidth', 2)
    ylim([-0.5 1])
    xticks(1:3:9);
    xticklabels(bar_name);
    for gl = 1:3
        disp([three_group_name{gl}, ':'])
        for bi = 1:3
            group_beta = sub_roi_beta(bi, rotation_label==gl)';
            [h,p,ci,stats] = ttest(group_beta);
            disp(['ttest, ', bar_name{bi}, '~=0: p-value=', num2str(p)])
            tvalue(gl, bi) = stats.tstat;
        end
        disp(' ')
    end




    %  Geometric relationships of VTC features across different contexts in the VTC on Day 3.
    roi_num = 5; % VTC
    if re_analysis
        cc_angle = (15:30:360)';
        cc_coord = [cosd(cc_angle), sind(cc_angle)];
        sub_feature_list = {};
        for sub_i = 1:size(sub_structure_model_objectEV, 1)
            sub_feature_list{sub_i, 1} = cc_coord;
        end
        [roi_features_dist, roi_rotation_rr] = s007_rdm_feature_rotation_and_compression(roi_cc_rdm, roi_num, sub_feature_list);
    else
        roi_cc_rotation_rr = load([figure_generation_data_path, 'r010_fMRI_ROIs_VTCfeatures_representational_geometry.mat']).roi_cc_rotation_rr;
        roi_rotation_rr = squeeze(roi_cc_rotation_rr(:, roi_num, :));
    end

    disp(' ')
    disp('figure : Geometric relationships of VTC features across different contexts in the VTC on Day 3')
    figure
    theta = deg2rad(1:361);
    group_color = {[255, 128, 0]./256, [7, 126, 151]./256};
    rotation_color = {[255, 216, 178]./256, [255, 164, 164]./256};
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
    rd_list = [-1, 1];
    rd_name = {'VTC_depended_clockwise', 'VTC_depended_anticlockwise'};
    for rd_i = 1:2
        rotation_rr = roi_rotation_rr(rotation_direction==rd_list(rd_i), :);
        SEM=std(rotation_rr)./sqrt(sum(rotation_direction==rd_list(rd_i)));
        [x1_SEM,y1_SEM] = pol2cart(theta, mean(rotation_rr)-SEM);
        [x2_SEM,y2_SEM] = pol2cart(theta, mean(rotation_rr)+SEM);
        fill([x1_SEM, fliplr(x2_SEM)], [y1_SEM, fliplr(y2_SEM)], 'k', 'FaceAlpha', ShadowAlpha, 'EdgeColor', 'none');
        hold on
        [x,y] = pol2cart(theta, mean(rotation_rr));
        plot(x, y, 'LineWidth', 1.7, 'color', rotation_color{rd_i})
        hold on
    end
    rotation_rr = roi_rotation_rr(group_label==2, :);
    SEM=std(rotation_rr)./sqrt(sum(group_label==2));
    [x1_SEM,y1_SEM] = pol2cart(theta, mean(rotation_rr)-SEM);
    [x2_SEM,y2_SEM] = pol2cart(theta, mean(rotation_rr)+SEM);
    fill([x1_SEM, fliplr(x2_SEM)], [y1_SEM, fliplr(y2_SEM)], 'k', 'FaceAlpha', ShadowAlpha, 'EdgeColor', 'none');
    hold on
    [x,y] = pol2cart(theta, mean(rotation_rr));
    plot(x, y, 'LineWidth', 1.7, 'color', group_color{2})
    hold on
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
    disp('Group level cross-context VTC features rotation in the VTC:')
    for gl = 1:2
        disp(group_name{gl})
        group_ra = sub_roi_ra(group_label==gl, 1);
        mean_rotation = mean(group_ra);
        SEM_rotation=std(group_ra)./sqrt(size(group_ra, 1));
        disp(['mean_rotation: ', num2str(mean_rotation)]);
        disp(['sem: ', num2str(SEM_rotation)])
        disp(' ')
    end
    for rd_i = 1:2
        disp(rd_name{rd_i})
        group_ra = sub_roi_ra(rotation_direction==rd_list(rd_i), 1);
        mean_rotation = mean(group_ra);
        SEM_rotation=std(group_ra)./sqrt(size(group_ra, 1));
        disp(['mean_rotation: ', num2str(mean_rotation)]);
        disp(['sem: ', num2str(SEM_rotation)])
        disp(' ')
    end



    % MDS
    % The figure may require rotation to match the results presented in the paper.
    disp(' ')
    disp('Visualization the Day 3 VTC representational geometry of object categories in a 2D plane through MDS.')
    original_rdm = [];
    roi_num = 5;
    for sub_i = 1:27
        roi_rdm = squareform(roi_cc_rdm{sub_i, roi_num, 4});
        original_rdm(:, :, sub_i) = roi_rdm;
    end
    angel = 15:30:360;
    center_colormap = hsv(360);
    center_color = center_colormap(angel, :);
    group_name = {'VTC_depended', 'Random_distribution'};
    context_color{1} = [0.4, 0.4, 0.4]; context_color{2} = [1, 0.15, 0.85]; context_color{3} = [0.15, 0.85, 1];
    for gi = 1:2
        disp(' ')
        disp(group_name{gi})
        group_average_rdm = squeeze(mean(original_rdm(:, :, group_label==gi), 3));
        [Y,e] = cmdscale(group_average_rdm, 2);
        
        figure
        disp('figure : indicate 12 categories with different color')
        scatter(Y(1:12, 1), Y(1:12, 2), 150, center_color, 'filled', "o"	)
        hold on
        scatter(Y(13:24, 1), Y(13:24, 2), 150, center_color, 'filled', "^")
        axis equal
        axis off

        figure
        disp('figure : indicate two context with different color')
        scatter(Y(1:12, 1), Y(1:12, 2), 150, context_color{2}, 'filled', "o"	)
        hold on
        scatter(Y(13:24, 1), Y(13:24, 2), 150, context_color{3}, 'filled', "^")
        axis equal
        axis off
    end
    rd_name = {'VTC_depended_clockwise', 'VTC_depended_anticlockwise'};
    for rd_i = 1:2
        disp(' ')
        disp(rd_name{rd_i})
        group_average_rdm = squeeze(mean(original_rdm(:, :, rotation_direction==rd_list(rd_i)), 3));
        [Y,e] = cmdscale(group_average_rdm, 2);

        figure
        disp('figure : indicate 12 categories with different color')
        scatter(Y(13:24, 1), Y(13:24, 2), 150, center_color, 'filled', "^")
        hold on
        scatter(Y(1:12, 1), Y(1:12, 2), 150, center_color, 'filled', "o"	)
        axis equal
        axis off
    end
    outputArg1 = 'Analysis completed without error';
end

