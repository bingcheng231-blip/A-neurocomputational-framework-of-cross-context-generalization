function [outputArg1] = f013_MEG_sensor_level_dynamic(Re_analysis_data_path,exp_label, figure_generation_data_path)
%F013_MEG_SENSOR_LEVEL_DYNAMIC 此处显示有关此函数的摘要
    [group_label, ~, ~, exp_index] = s002_grouplabel(Re_analysis_data_path, exp_label);
 
    sub_localizer_tp_rdm = load([Re_analysis_data_path, 'd014_MEG_sensor_space_tp_rdm_localizer.mat']).sub_localizer_tp_rdm;
    sub_tp_rdm = [];
    for sub_i = 1:23
        sub_tp_rdm(sub_i, :, :) = sub_localizer_tp_rdm{sub_i, 1};
    end

    sub_corr = [];
    for sub_i = 1:23
        others_lb = setdiff(1:23, sub_i);
        sub_mean = squeeze(mean(sub_tp_rdm(others_lb, :, :)));
        for ti = 1:151
            sub_corr(sub_i, ti) = corr(squeeze(sub_tp_rdm(sub_i, :, ti))', squeeze(sub_mean(:, ti)));
        end
    end
    % stimulus stage
    disp(' ')
    disp('figure : MEG localizer noise ceiling')
    MEG_noise_ceiling = sub_corr(:, 34:end-17);
    time = linspace(-200, 800, 101);
    group_name = {'VTC_dependent', 'Random_distribution'};
    for gl = 1:2
        disp(group_name{gl})
        figure
        x=time;
        y = squeeze(MEG_noise_ceiling(group_label==gl, :));
        SEM=std(y)./sqrt(size(y, 1));
        s = shadedErrorBar(x,mean(y, 1),SEM, 'lineProps', '-');
        s.patch.FaceColor = [0.5, 0.5, 0.5];
        set(s.edge,'LineWidth',0.5, 'Color', [0.7, 0.7, 0.7])
        hold on
        plot(time, mean(y, 1), 'LineWidth', 1.5, 'color', ones(1, 3).*0.5)
    end

    cc_angle = (15:30:360)';
    cc_coord = [cosd(cc_angle), sind(cc_angle)];
    cc_rdm = zscore(pdist(cc_coord))';
    model_beta = [];
    for sub_i = 1:23
        for ti = 1:151
            model_beta(sub_i, ti) = corr(squeeze(sub_tp_rdm(sub_i, :, ti))', cc_rdm);
        end
    end
    model_beta = model_beta(:, 34:end-17);
    disp(' ')
    disp('figure : category representation during MEG localizer')
    for gl = 1:2
        disp(group_name{gl})
        figure
        x=time;
        y = squeeze(model_beta(group_label==1, :));
        SEM=std(y)./sqrt(size(y, 1));
        s = shadedErrorBar(x,mean(y, 1),SEM, 'lineProps', '-');
        s.patch.FaceColor = [0.5, 0.5, 0.5];
        set(s.edge,'LineWidth',0.5, 'Color', [0.7, 0.7, 0.7])
        hold on
        plot(time, mean(y, 1), 'LineWidth', 1.5, 'color', ones(1, 3).*0.5)
    end



    fMRI_sub_NO = nan(27, 1);
    sub_info = load([Re_analysis_data_path, 'd004_Sub_learning_trials_info.mat']).sub_learning_trials_info;
    for sub_i = 1:36
        exp_index(sub_i, 1) = contains(sub_info{sub_i, 1}, 'fMRI');
    end
    fMRI_sub_list = sub_info(exp_index, 1);
    for sub_i = 1:27
        fMRI_sub_NO(sub_i) = str2num(fMRI_sub_list{sub_i}(8:9));
    end

    MEG_sub_NO = nan(23, 1);
    sub_info = load([Re_analysis_data_path, 'd004_Sub_learning_trials_info.mat']).sub_learning_trials_info;
    for sub_i = 1:36
        exp_index(sub_i, 1) = contains(sub_info{sub_i, 1}, 'MEG');
    end
    MEG_sub_list = sub_info(exp_index, 1);
    for sub_i = 1:23
        MEG_sub_NO(sub_i) = str2num(MEG_sub_list{sub_i}(8:9));
    end
    %  Geometric relationships of EVs across different contexts in ROIs on Day 3.
    STRUCT_flag = true;
    all_behavior_model = load([Re_analysis_data_path, 'd007_behavior_model.mat']).all_behavior_model;
    [sub_structure_model_objectEV] = s005_expected_value_of_object(Re_analysis_data_path, 'fMRI', STRUCT_flag, all_behavior_model);

    cc_geometry = load([figure_generation_data_path, 'r010_fMRI_ROIs_VTCfeatures_representational_geometry.mat']);
    vv_geometry = load([figure_generation_data_path, 'r009_fMRI_ROIs_EV_representational_geometry.mat']);
    ROI_num = 2;
    vv_rdm = []; cc_rdm = []; vv_param = {};
    for sub_i = 1:27
        structure_EVs = sub_structure_model_objectEV{sub_i, 1}(end, :);
        sub_vv_coord = [structure_EVs(1, 1:12)', structure_EVs(1, 13:24)'];
        [M,I] = max(squeeze(vv_geometry.roi_vv_rotation_rr(sub_i, ROI_num, :)));
        rotation = [[cosd(I); sind(I)], [-sind(I); cosd(I)]];
        vv_param{sub_i, 1} = rotation;
        vv_param{sub_i, 2} = vv_geometry.roi_ctx_dist(sub_i, ROI_num, I, :);
        vmPFC_vv_rdm = s008_ROIs_rdm_ctxDist_calc(sub_vv_coord, rotation, vv_geometry.roi_ctx_dist(sub_i, ROI_num, I, :));
        vv_rdm(:, sub_i) = zscore(vmPFC_vv_rdm');

        [M,I] = max(squeeze(cc_geometry.roi_cc_rotation_rr(sub_i, ROI_num, :)));
        rotation = [[cosd(I); sind(I)], [-sind(I); cosd(I)]];
        vmPFC_cc_rdm = s008_ROIs_rdm_ctxDist_calc(cc_coord, rotation, cc_geometry.roi_vtc_features_dist(sub_i, ROI_num, I, :));
        cc_rdm(:, sub_i) = zscore(vmPFC_cc_rdm');
    end

    MEG_sub_learning_stimulus_tp_rdm = load(fullfile(Re_analysis_data_path, 'd015_MEG_sub_learning_stimulus_tp_rdm.mat')).MEG_sub_learning_stimulus_tp_rdm;
    sub_roi_beta = [];
    for sub_i = 1:23
        sub_NO = MEG_sub_NO(sub_i);

        sub_vv_rdm = vv_rdm(:, fMRI_sub_NO==sub_NO);
        if sum(fMRI_sub_NO==sub_NO) > 1
            sub_vv_rdm = mean(sub_vv_rdm, 2);
        end
        sub_vv_rdm = zscore(sub_vv_rdm);

        sub_cc_rdm = cc_rdm(:, fMRI_sub_NO==sub_NO);
        if sum(fMRI_sub_NO==sub_NO) > 1
            sub_cc_rdm = mean(sub_cc_rdm, 2);
        end
        sub_cc_rdm = zscore(sub_cc_rdm);

        ctx_dist = [zeros(12, 1); ones(12, 1)];
        structure_EVs = sub_structure_model_objectEV{sub_i, 1}(end, :);
        sub_vv_coord = [structure_EVs(1, 1:12)', structure_EVs(1, 13:24)'];
        sub_vv_rdm(:, 2)  = zscore(pdist([[sub_vv_coord; sub_vv_coord], ctx_dist]));
        sub_vv_rdm(:, 3)  = zscore(pdist([[[sub_vv_coord(:, 1), zeros(12, 1)]; [zeros(12, 1), sub_vv_coord(:, 2)]], ctx_dist]));
        sub_vv_rdm(:, 4)  = zscore(pdist([[[sub_vv_coord(:, 1), zeros(12, 1)]; [sub_vv_coord(:, 2), zeros(12, 1)]], ctx_dist]));

        sub_cc_rdm(:, 2)  = zscore(pdist([[cc_coord; cc_coord], ctx_dist]));
        sub_cc_rdm(:, 3)  = zscore(pdist([[[cc_coord(:, 1), zeros(12, 1)]; [zeros(12, 1), cc_coord(:, 2)]], ctx_dist]));
        sub_cc_rdm(:, 4)  = zscore(pdist([[[cc_coord(:, 1), zeros(12, 1)]; [cc_coord(:, 2), zeros(12, 1)]], ctx_dist]));

        model_rdm = [];
        model_rdm = [sub_vv_rdm, sub_cc_rdm];

        baseline_rdm = zscore(MEG_sub_learning_stimulus_tp_rdm{sub_i, 2});
        target_roi_rdm = zscore(MEG_sub_learning_stimulus_tp_rdm{sub_i, 1});
        for ti = 1:131
            mdl_beta = lsqnonneg([model_rdm, baseline_rdm, ones(276, 1)], target_roi_rdm(:, ti));
            sub_roi_beta(sub_i, :, ti) = mdl_beta(1:end-1);
        end
    end
    % stimulus stage
    sub_vc_beta{2} = sub_roi_beta(:, [5:8, 9], 34:end-17); % from triggel to MEG visual need 33ms
    sub_vc_beta{1} = sub_roi_beta(:, [1:4, 9], 34:end-17); % from triggel to MEG visual need 33ms
    time = linspace(-200, 600, 81);

    vc_color{1, 1} = [0.9290 0.5940 0.1250]; % value
    vc_color{1, 2} = [0.5000 0.0000 0.2250]; % value
    vc_color{1, 3} = [1.0000 0.9000 0.3250]; % value
    vc_color{1, 4} = [1.0000 0.3333 0.0000]; % value

    vc_color{2, 1} = [0  , 156, 0 ]./255; % category
    vc_color{2, 2} = [0  , 255, 64 ]./255; % category
    vc_color{2, 3} = [0  , 156, 156 ]./255; % category
    vc_color{2, 4} = [0  , 0.9333,   1.0000]; % category

    disp(' ')
    disp('figure : Temporal Evolution of MEG Sensor-Level Activity During Stimulus Presentation')
    group_name = {'VTC_dependent', 'Random_distribution'};
    vc_name = {'value info', 'category info'};
    for vc_num = 1:2
        disp(' ')
        disp(vc_name{vc_num})
        cluster_forming_threshold = 0.05;
        rps = RandStream('dsfmt19937');
        for gl = 1:2
            disp(group_name{gl})
            figure
            for beta_i = 1:4
                y = squeeze(sub_vc_beta{vc_num}(group_label==gl, beta_i, :));
                b = squeeze(sub_vc_beta{vc_num}(group_label==gl, end, :));

                SEM=std(y)./sqrt(size(y, 1));
                s = shadedErrorBar(time,mean(y, 1),SEM, 'lineProps', '-');
                s.patch.FaceColor = [0.5, 0.5, 0.5];
                set(s.edge,'LineWidth',0.5, 'Color', [0.7, 0.7, 0.7])
                hold on
                plot(time, mean(y, 1), 'LineWidth', 1.5, 'color',  vc_color{vc_num, beta_i})

                [h,p,ci,stats] = ttest(y, b, 'tail', 'right');
                tp_pvalue = p;
                tp_bw = tp_pvalue<cluster_forming_threshold;
                tp_CC = bwconncomp(tp_bw);
                tp_cluster_size = nan(tp_CC.NumObjects, 1);
                for ci = 1:tp_CC.NumObjects
                    tp_cluster_size(ci) = numel(tp_CC.PixelIdxList{ci});
                end
                tp_cluster = tp_CC.PixelIdxList;

                % cluster correction
                perm_num = 10000;
                cluster_size = nan(perm_num, 1);
                all_data = [y; b]; half_n = size(all_data, 1)./2;
                parfor perm_i = 1:perm_num
                    sub_prem = randperm(rps, size(all_data, 1), size(all_data, 1))';
                    prem_data = all_data(sub_prem, :);
                    perm_y = prem_data(1:half_n, :);
                    perm_b = prem_data(half_n+1:end, :);
                    [~,perm_p,~,~] = ttest(perm_y, perm_b, 'tail', 'right');

                    perm_bw = perm_p<cluster_forming_threshold;
                    perm_CC = bwconncomp(perm_bw);
                    perm_cluster = nan(perm_CC.NumObjects, 1);
                    if perm_CC.NumObjects>0
                        for ci = 1:perm_CC.NumObjects
                            perm_cluster(ci) = numel(perm_CC.PixelIdxList{ci});
                        end
                        max_size = max(perm_cluster);
                        cluster_size(perm_i) = max_size(1);
                    else
                        cluster_size(perm_i) = 0;
                    end
                end
                sorted_cluster_sizes = sort(cluster_size(~isnan(cluster_size)),'descend');

                pthresh=0.05;
                cluster_size_thresh = round(interp1(log2((1:perm_num)/perm_num), sorted_cluster_sizes, log2(pthresh)));
                unsurvival_tp_cluster = tp_cluster(tp_cluster_size<cluster_size_thresh);
                for sc_i = 1:length(unsurvival_tp_cluster)
                    tp_pvalue(unsurvival_tp_cluster{sc_i}) = 1;
                end

                p_line = nan(size(tp_pvalue));
                p_line(tp_pvalue<cluster_forming_threshold) = -0.05-0.01.*beta_i;
                plot(time, p_line, 'LineWidth', 2, 'color',  vc_color{vc_num, beta_i})
                plot(time, mean(y, 1), 'LineWidth', 1.5, 'color',  vc_color{vc_num, beta_i})
            end
            hold on
            plot([min(time), max(time)],[0, 0], '--k', 'LineWidth', 1)
            hold on
            plot([0, 0], [-0.2, 0.6], '--k', 'LineWidth', 1)
            set(gca, 'FontSize', 15, 'FontWeight', 'bold')
            set(gca,'LineWidth',2)
            set(gca, 'position', [0.12, 0.15, 0.8, 0.7]);
        end
    end



    outputArg1 = 'Analysis completed without error';
end

