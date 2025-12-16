function [outputArg1] = f014_MEG_ROI_interactions(Re_analysis_data_path,exp_label, figure_generation_data_path)
%F014_MEG_ROI_INTERACTIONS 此处显示有关此函数的摘要
%   此处显示详细说明
    it4=wd_compit4_1;

    disp(' ')
    disp('figure : ROI dynamic interactions during the MEG localizer experiment  ')
    disp('It takes some time to generate the result images....')
    prem_cluster_correction = true;
    clip_time = 114; % stimulus stage
    began_tp = 34; % stimulus stage
    past_deltaT = 0;
    roi_tp_rdm = load([Re_analysis_data_path, 'd016_MEG_ROIs_source_space_tp_rdm_localizer.mat']).roi_tp_rdm;
    roi_name = {'dlPFC', 'vmPFC', 'OFC', 'MCC', 'EC', 'VTC'};
    for source_ROInum = 2:6
        for target_ROInum = setdiff(2:6, source_ROInum)
            sub_beta = zeros(clip_time, clip_time, 2, 23);
            for sub_i = 1:23
                target_roi_rdm = roi_tp_rdm{sub_i, target_ROInum};
                target_roi_rdm = zscore(target_roi_rdm);

                for past_ti = began_tp:clip_time
                    target_past = zscore(mean(roi_tp_rdm{sub_i, target_ROInum}(:, past_ti-past_deltaT:past_ti), 2));
                    target_base = zscore(mean(roi_tp_rdm{sub_i, target_ROInum}(:, 34:53), 2));
                    target_model = [target_past, target_base];

                    source_past = zscore(mean(roi_tp_rdm{sub_i, source_ROInum}(:, past_ti-past_deltaT:past_ti), 2));
                    source_base = zscore(mean(roi_tp_rdm{sub_i, source_ROInum}(:, 34:53), 2));
                    for tp_i = began_tp:clip_time
                        mdl = lsqnonneg([source_past, source_base, target_model, ones(66, 1)], target_roi_rdm(:, tp_i));
                        sub_beta(clip_time+1-past_ti, tp_i, 1:2, sub_i) = mdl(1:2);
                    end
                end
            end

            south_triangle = tril(true(clip_time), -1);
            south_triangle = rot90(south_triangle);
            south_triangle(clip_time-began_tp+2:clip_time, :) = false;
            south_beta = 0.05.*ones(clip_time, clip_time, 23);
            baseline_beta = 0.05.*ones(clip_time, clip_time, 23);
            line_beta = []; baseline = [];
            for sub_i = 1:23
                tp_beta = sub_beta(:, :, 1, sub_i);
                line_beta(sub_i, :) = tp_beta(south_triangle);
                tp_beta(~south_triangle)=0.05;
                south_beta(:, :, sub_i) = tp_beta;

                tp_beta = sub_beta(:, :, 2, sub_i);
                baseline(sub_i, :) = tp_beta(south_triangle);
                baseline_beta(:, :, sub_i) = tp_beta;
            end
            south_aps = mean(south_beta-baseline_beta, 3);

            [h,p,ci,stats] = ttest(line_beta, baseline, 'tail', 'right');
            tp_pvalue = ones(clip_time, clip_time);
            tp_pvalue(south_triangle) = p;
            tp_bw = tp_pvalue<0.05;
            tp_CC = bwconncomp(tp_bw);
            tp_cluster_size = nan(tp_CC.NumObjects, 1);
            for ci = 1:tp_CC.NumObjects
                tp_cluster_size(ci) = numel(tp_CC.PixelIdxList{ci});
            end
            tp_cluster = tp_CC.PixelIdxList;

            if prem_cluster_correction
                % cluster_correction
                perm_num = 10000;
                cluster_size = nan(perm_num, 1);
                all_data = [line_beta; baseline];
                parfor perm_i = 1:perm_num
                    sub_prem = randperm(23*2, 23*2)';
                    perm_data = all_data(sub_prem, :);
                    [~,perm_p,~,~] = ttest(perm_data(1:23, :), perm_data(24:end, :), 'tail', 'right');
                    % [~,perm_p,~,~] = ttest(perm_data(1:23, :), perm_data(24:end, :));
                    perm_pvalue = ones(clip_time, clip_time);
                    perm_pvalue(south_triangle) = perm_p;
                    perm_bw = perm_pvalue<0.05;
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
                % pthresh=0.01./comp_right.NumObjects; % FWE-corrected P values
                cluster_size_thresh = round(interp1(log2((1:perm_num)/perm_num), sorted_cluster_sizes, log2(pthresh)));
                unsurvival_tp_cluster = tp_cluster(tp_cluster_size<cluster_size_thresh);
                for sc_i = 1:length(unsurvival_tp_cluster)
                    tp_pvalue(unsurvival_tp_cluster{sc_i}) = 1;
                end
            end
            south_tp_pvalue = tp_pvalue;

            sub_beta = zeros(clip_time, clip_time, 2, 23);
            for sub_i = 1:23
                target_roi_rdm = roi_tp_rdm{sub_i, source_ROInum};
                target_roi_rdm = zscore(target_roi_rdm);
                for past_ti = began_tp:clip_time
                    target_past = zscore(mean(roi_tp_rdm{sub_i, source_ROInum}(:, past_ti-past_deltaT:past_ti), 2));
                    target_base = zscore(mean(roi_tp_rdm{sub_i, source_ROInum}(:, 34:53), 2));
                    target_model = [target_past, target_base];
                    source_past = zscore(mean(roi_tp_rdm{sub_i, target_ROInum}(:, past_ti-past_deltaT:past_ti), 2));
                    source_base = zscore(mean(roi_tp_rdm{sub_i, target_ROInum}(:, 34:53), 2));
                    for tp_i = began_tp:clip_time
                        mdl = lsqnonneg([source_past, source_base, target_model, ones(66, 1)], target_roi_rdm(:, tp_i));
                        sub_beta(clip_time+1-past_ti, tp_i, 1:2, sub_i) = mdl(1:2);
                    end
                end
            end
            north_triangle = triu(true(clip_time), 1);
            north_triangle = rot90(north_triangle);
            north_triangle(:, 1:began_tp-1) = false;
            north_beta = 0.05.*ones(clip_time, clip_time, 23);
            baseline_beta = 0.05.*ones(clip_time, clip_time, 23);
            line_beta = []; baseline = [];
            for sub_i = 1:23
                tp_beta = sub_beta(:, :, 1, sub_i);
                line_beta(sub_i, :) = tp_beta(north_triangle);
                tp_beta(~north_triangle)=0.05;
                north_beta(:, :, sub_i) = tp_beta;

                tp_beta = sub_beta(:, :, 2, sub_i);
                baseline(sub_i, :) = tp_beta(south_triangle);
                baseline_beta(:, :, sub_i) = tp_beta;
            end
            north_aps = mean(north_beta-baseline_beta, 3);

            [h,p,ci,stats] = ttest(line_beta, baseline, 'tail', 'right');
            tp_pvalue = ones(clip_time, clip_time);
            tp_pvalue(south_triangle) = p;
            tp_bw = tp_pvalue<0.05;
            tp_CC = bwconncomp(tp_bw);
            tp_cluster_size = nan(tp_CC.NumObjects, 1);
            for ci = 1:tp_CC.NumObjects
                tp_cluster_size(ci) = numel(tp_CC.PixelIdxList{ci});
            end
            tp_cluster = tp_CC.PixelIdxList;

            if prem_cluster_correction
                % cluster_correction
                perm_num = 10000;
                cluster_size = nan(perm_num, 1);
                all_data = [line_beta; baseline];
                parfor perm_i = 1:perm_num
                    sub_prem = randperm(23*2, 23*2)';
                    perm_data = all_data(sub_prem, :);
                    [~,perm_p,~,~] = ttest(perm_data(1:23, :), perm_data(24:end, :), 'tail', 'right');
                    % [~,perm_p,~,~] = ttest(perm_data(1:23, :), perm_data(24:end, :));
                    perm_pvalue = ones(clip_time, clip_time);
                    perm_pvalue(south_triangle) = perm_p;
                    perm_bw = perm_pvalue<0.05;
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
                % pthresh=0.01./comp_right.NumObjects; % FWE-corrected P values
                cluster_size_thresh = round(interp1(log2((1:perm_num)/perm_num), sorted_cluster_sizes, log2(pthresh)));
                unsurvival_tp_cluster = tp_cluster(tp_cluster_size<cluster_size_thresh);
                for sc_i = 1:length(unsurvival_tp_cluster)
                    tp_pvalue(unsurvival_tp_cluster{sc_i}) = 1;
                end
            end
            north_tp_pvalue = tp_pvalue;

            aps = 0.5.*ones(clip_time, clip_time);
            aps(south_triangle) = (south_aps(south_triangle)+0.1)./0.2;
            aps(north_triangle) = (north_aps(north_triangle)+0.1)./0.2;
            for ti = 4:10:clip_time
                aps(ti-3, began_tp-1) = 0;
                aps(clip_time-began_tp+2, ti) = 0;
            end
            map_corr=[];
            for i=1:size(aps,1)
                for j=1:size(aps,2)
                    map_corr(i,j,1:3)=it4(max(min(round(255*aps(i,j))+1,256),1),1:3);
                end
            end
            f = figure;
            f.Position(3:4) = [500 500];
            imshow(map_corr, 'border', 'tight', 'InitialMagnification','fit');
            hold on

            tp_pvalue = ones(clip_time, clip_time);
            tp_pvalue(south_triangle) = south_tp_pvalue(south_triangle);
            tp_pvalue(north_triangle) = north_tp_pvalue(north_triangle);
            [X, Y] = meshgrid(1:clip_time);
            v = [0.05, 0.05];
            contour(X, Y, tp_pvalue, v, 'LineWidth', 1.5, 'LineColor', 'g')
            hold on
            plot([began_tp+20 began_tp+20], [1 clip_time-began_tp+2], '-k', 'LineWidth', 1.5)
            plot([began_tp-1 clip_time],[clip_time-began_tp-19 clip_time-began_tp-19],  '-k', 'LineWidth', 1.5)
            plot([clip_time began_tp],[1 clip_time-began_tp+1],  '-k', 'LineWidth', 1.5)
        end
    end



    disp(' ')
    disp('figure : ROI dynamic interactions during the MEG learning experiment  ')
    disp('It takes some time to generate the result images....')
    clip_time = 114; % stimulus stage
    began_tp = 34; % stimulus stage
    past_deltaT = 0;
    cluster_forming = 0.05;
    pthresh=0.05;
    prem_cluster_correction = true;

    sub_roi_tp_rdm = load(fullfile(Re_analysis_data_path, 'd017_MEG_ROIs_source_space_tp_rdm_learning_stimulus.mat')).sub_roi_tp_rdm;
    sub_roi_base_line = load(fullfile(Re_analysis_data_path, 'd017_MEG_ROIs_source_space_tp_rdm_learning_stimulus.mat')).sub_roi_base_line;
    roi_name = {'dlPFC', 'vmPFC', 'OFC', 'MCC', 'EC', 'VTC'};
    rps = RandStream('dsfmt19937');
    for source_ROInum = 2:6
        for target_ROInum = setdiff(2:6, source_ROInum)
            sub_beta = zeros(clip_time, clip_time, 2, 23);
            for sub_i = 1:23
                base_line = sub_roi_base_line(:, sub_i);
                roi_tp_rdm = sub_roi_tp_rdm(:, sub_i);

                target_roi_rdm = roi_tp_rdm{target_ROInum, 1};
                target_roi_rdm = zscore(target_roi_rdm);

                for past_ti = began_tp:clip_time
                    target_past = zscore(mean(roi_tp_rdm{target_ROInum, 1}(:, past_ti-past_deltaT:past_ti), 2));
                    target_base = zscore(mean(base_line{target_ROInum, 1}(:, 14:33), 2));
                    target_model = [target_past, target_base];

                    source_past = zscore(mean(roi_tp_rdm{source_ROInum, 1}(:, past_ti-past_deltaT:past_ti), 2));
                    source_base = zscore(mean(base_line{source_ROInum, 1}(:, 14:33), 2));
                    for tp_i = began_tp:clip_time
                        mdl = lsqnonneg([source_past, source_base, target_model, ones(276, 1)], target_roi_rdm(:, tp_i));
                        sub_beta(clip_time+1-past_ti, tp_i, 1:2, sub_i) = mdl(1:2);
                    end
                end
            end
            south_triangle = tril(true(clip_time), -1);
            south_triangle = rot90(south_triangle);
            south_triangle(clip_time-began_tp+2:clip_time, :) = false;
            south_beta = 0.05.*ones(clip_time, clip_time, size(sub_beta, 4));
            south_base = 0.05.*ones(clip_time, clip_time, size(sub_beta, 4));
            line_beta = []; baseline = [];
            for sub_i = 1:size(sub_beta, 4)
                tp_beta = sub_beta(:, :, 1, sub_i);
                line_beta(sub_i, :) = tp_beta(south_triangle);
                tp_beta(~south_triangle)=0.05;
                south_beta(:, :, sub_i) = tp_beta;

                tp_beta = sub_beta(:, :, 2, sub_i);
                baseline(sub_i, :) = tp_beta(south_triangle);
                south_base(:, :, sub_i) = tp_beta;
            end
            south_aps = mean(south_beta-south_base, 3);

            [h,p,ci,stats] = ttest(line_beta, baseline, 'tail', 'right');
            tp_pvalue = ones(clip_time, clip_time);
            tp_pvalue(south_triangle) = p;
            tp_bw = tp_pvalue<cluster_forming;
            tp_CC = bwconncomp(tp_bw);
            tp_cluster_size = nan(tp_CC.NumObjects, 1);
            for ci = 1:tp_CC.NumObjects
                tp_cluster_size(ci) = numel(tp_CC.PixelIdxList{ci});
            end
            tp_cluster = tp_CC.PixelIdxList;

            if prem_cluster_correction
                % cluster_correction
                perm_num = 10000;
                cluster_size = nan(perm_num, 1);
                all_data = [line_beta; baseline];
                parfor perm_i = 1:perm_num
                    sub_prem = randperm(rps, size(sub_beta, 4)*2, size(sub_beta, 4)*2)';
                    perm_data = all_data(sub_prem, :);
                    [~,perm_p,~,~] = ttest(perm_data(1:size(sub_beta, 4), :), perm_data(size(sub_beta, 4)+1:end, :), 'tail', 'right');
                    perm_pvalue = ones(clip_time, clip_time);
                    perm_pvalue(south_triangle) = perm_p;
                    perm_bw = perm_pvalue<cluster_forming;
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
                cluster_size_thresh = round(interp1(log2((1:perm_num)/perm_num), sorted_cluster_sizes, log2(pthresh)));
                unsurvival_tp_cluster = tp_cluster(tp_cluster_size<cluster_size_thresh);
                for sc_i = 1:length(unsurvival_tp_cluster)
                    tp_pvalue(unsurvival_tp_cluster{sc_i}) = 1;
                end
            end
            south_tp_pvalue = tp_pvalue;
            south_tp_tvalue = ones(clip_time, clip_time);
            south_tp_tvalue(south_triangle) = stats.tstat;

            sub_beta = zeros(clip_time, clip_time, 2, 23);
            for sub_i = 1:23
                base_line = sub_roi_base_line(:, sub_i);
                roi_tp_rdm = sub_roi_tp_rdm(:, sub_i);

                target_roi_rdm = roi_tp_rdm{source_ROInum, 1};
                target_roi_rdm = zscore(target_roi_rdm);
                for past_ti = began_tp:clip_time
                    target_past = zscore(mean(roi_tp_rdm{source_ROInum, 1}(:, past_ti-past_deltaT:past_ti), 2));
                    target_base = zscore(mean(base_line{source_ROInum, 1}(:, 14:33), 2));
                    target_model = [target_past, target_base];

                    source_past = zscore(mean(roi_tp_rdm{target_ROInum, 1}(:, past_ti-past_deltaT:past_ti), 2));
                    source_base = zscore(mean(base_line{target_ROInum, 1}(:, 14:33), 2));
                    for tp_i = began_tp:clip_time
                        mdl = lsqnonneg([source_past, source_base, target_model, ones(276, 1)], target_roi_rdm(:, tp_i));
                        sub_beta(clip_time-tp_i+1, past_ti, 1:2, sub_i) = mdl(1:2);
                    end
                end
            end
            north_triangle = triu(true(clip_time), 1);
            north_triangle = rot90(north_triangle);
            north_triangle(:, 1:began_tp-1) = false;
            north_beta = 0.05.*ones(clip_time, clip_time, size(sub_beta, 4));
            north_base = 0.05.*ones(clip_time, clip_time, size(sub_beta, 4));
            line_beta = []; baseline = [];
            for sub_i = 1:size(sub_beta, 4)
                tp_beta = sub_beta(:, :, 1, sub_i);
                line_beta(sub_i, :) = tp_beta(north_triangle);
                tp_beta(~north_triangle)=0.05;
                north_beta(:, :, sub_i) = tp_beta;

                tp_beta = sub_beta(:, :, 2, sub_i);
                baseline(sub_i, :) = tp_beta(north_triangle);
                north_base(:, :, sub_i) = tp_beta;
            end
            north_aps = mean(north_beta-north_base, 3);

            [h,p,ci,stats] = ttest(line_beta, baseline, 'tail', 'right');
            tp_pvalue = ones(clip_time, clip_time);
            tp_pvalue(north_triangle) = p;
            tp_bw = tp_pvalue<cluster_forming;
            tp_CC = bwconncomp(tp_bw);
            tp_cluster_size = nan(tp_CC.NumObjects, 1);
            for ci = 1:tp_CC.NumObjects
                tp_cluster_size(ci) = numel(tp_CC.PixelIdxList{ci});
            end
            tp_cluster = tp_CC.PixelIdxList;

            if prem_cluster_correction
                % cluster_correction
                perm_num = 10000;
                cluster_size = nan(perm_num, 1);
                all_data = [line_beta; baseline];
                parfor perm_i = 1:perm_num
                    sub_prem = randperm(rps, size(sub_beta, 4)*2, size(sub_beta, 4)*2)';
                    perm_data = all_data(sub_prem, :);
                    [~,perm_p,~,~] = ttest(perm_data(1:size(sub_beta, 4), :), perm_data(size(sub_beta, 4)+1:end, :), 'tail', 'right');
                    perm_pvalue = ones(clip_time, clip_time);
                    perm_pvalue(north_triangle) = perm_p;
                    perm_bw = perm_pvalue<cluster_forming;
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
                cluster_size_thresh = round(interp1(log2((1:perm_num)/perm_num), sorted_cluster_sizes, log2(pthresh)));
                unsurvival_tp_cluster = tp_cluster(tp_cluster_size<cluster_size_thresh);
                for sc_i = 1:length(unsurvival_tp_cluster)
                    tp_pvalue(unsurvival_tp_cluster{sc_i}) = 1;
                end
            end
            north_tp_pvalue = tp_pvalue;
            north_tp_tvalue = ones(clip_time, clip_time);
            north_tp_tvalue(north_triangle) = stats.tstat;

            aps = 0.5.*ones(clip_time, clip_time);
            aps(south_triangle) = (south_aps(south_triangle)+0.1)./0.2;
            aps(north_triangle) = (north_aps(north_triangle)+0.1)./0.2;
            for ti = 4:10:clip_time
                aps(clip_time-ti+1, began_tp-1) = 0;
                aps(clip_time-began_tp+2, ti) = 0;
            end
            map_corr=[];
            for i=1:size(aps,1)
                for j=1:size(aps,2)
                    map_corr(i,j,1:3)=it4(max(min(round(255*aps(i,j))+1,256),1),1:3);
                end
            end
            f = figure;
            f.Position(3:4) = [500 500];
            imshow(map_corr, 'border', 'tight', 'InitialMagnification','fit');
            hold on

            tp_pvalue = ones(clip_time, clip_time);
            tp_pvalue(south_triangle) = south_tp_pvalue(south_triangle);
            tp_pvalue(north_triangle) = north_tp_pvalue(north_triangle);
            [X, Y] = meshgrid(1:clip_time);
            v = [cluster_forming, cluster_forming];
            contour(X, Y, tp_pvalue, v, 'LineWidth', 1.5, 'LineColor', 'g')
            hold on
            plot([began_tp+20 began_tp+20], [1 clip_time-began_tp+2], '-k', 'LineWidth', 1.5)
            plot([began_tp-1 clip_time],[clip_time-began_tp-19 clip_time-began_tp-19],  '-k', 'LineWidth', 1.5)
            plot([clip_time began_tp],[1 clip_time-began_tp+1],  '-k', 'LineWidth', 1.5)
        end
    end
    outputArg1 = 'Analysis completed without error';
end

