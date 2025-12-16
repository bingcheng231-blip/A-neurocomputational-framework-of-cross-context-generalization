function [outputArg1] = f006_vmPFC_model_RDM_fitting(Re_analysis_data_path,exp_label)
%F006_VMPFC_MODEL_RDM_FITTING 此处显示有关此函数的摘要
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
    ctx_vv_beta = nan(27, 3, 3);
    roi_num = 2; % vmPFC
    for sub_i = 1:27
        structure_model_EV(1, :) = sub_structure_model_objectEV{sub_i, 1}(end, :);
        sub_vv_coord = [structure_model_EV(1, 1:12)', structure_model_EV(1, 13:24)'];
        ctxA_vv_rdm = zscore(pdist(sub_vv_coord(:, 1)))';
        ctxB_vv_rdm = zscore(pdist(sub_vv_coord(:, 2)))';
        for ctx_i = 1:3
            roi_rdm = zscore(roi_cc_rdm{sub_i, roi_num, ctx_i})';
            if ~isempty(roi_rdm)
                ctx_vv_beta(sub_i, ctx_i, :) = regress(roi_rdm, [ctxA_vv_rdm, ctxB_vv_rdm, ones(size(ctxB_vv_rdm))]);
            end
        end
    end


    group_vv_beta = cell(2, 3);
    for gl = 1:2
        for ctx_i = 1:3
            vv_beta = squeeze(ctx_vv_beta(group_label==gl, ctx_i, 1:end-1));
            group_vv_beta{gl, ctx_i} = vv_beta;
        end
    end

    disp(['figure : Representations of EVs in the vmPFC across different days and contexts'])
    disp(' ')
    figure;
    group_color = {[255, 128, 0]./256, [7, 126, 151]./256};
    group_name = {'VTC_depended_group', 'Random_distribution_group'};
    bn = 1;
    bar_name = {'day01 EV01', 'day01 EV02', 'day03 CtxA EV01', 'day03 CtxA EV02', ...
         'day03 CtxB EV01', 'day03 CtxB EV02'};
    for ctx_i = 1:3
        for fi = 1:2
            for gl = 1:2
                gm_beta = group_vv_beta{gl, ctx_i}(:, fi);
                gm_beta = gm_beta(~isnan(gm_beta));
                gm_means = mean(gm_beta);
                gm_sems  = std(gm_beta)/sqrt(size(gm_beta, 1));
                jitter = 0.2;
                b = bar(bn-0.5, gm_means, 0.6, 'FaceColor','none', 'EdgeColor','none');
                hold on
                % xi = bn-0.5 + (rand(size(gm_beta)) - 0.5) * 2 * jitter;
                % scatter(xi, gm_beta, 50, ones(1, 3).*0.2, 'filled', ...
                %     'MarkerFaceAlpha',0.7, 'MarkerEdgeAlpha',0.2);
                b.EdgeColor = group_color{gl};
                b.LineWidth = 5;
                errorbar(bn-0.5, gm_means, gm_sems, 'color', group_color{gl}, 'LineWidth',2, 'CapSize',12, 'LineStyle','none');
                bn = bn+1;
            end
        end
    end
    box off
    set(gca, 'FontSize', 15, 'FontWeight', 'bold')
    set(gca,'LineWidth', 2)
    ylim([-0.2 0.8])
    xticks(1:2:14);
    xticklabels(bar_name);

    for gl = 1:2
        disp([group_name{gl}, ':'])
        group_beta = cell2mat(group_vv_beta(gl, :));
        for bi = 1:size(group_beta, 2)
            [h,p,ci,stats] = ttest(group_beta(:, bi), 0, 'tail', 'right');
            disp(['ttest, ', bar_name{bi}, '>0: p-value=', num2str(p)])
            tvalue(gl, bi) = stats.tstat;
        end
        disp(' ')
    end

  
    % EV model comparison through RSA
    sub_roi_beta = [];
    for sub_i = 1:27
        structure_model_EV(1, :) = sub_structure_model_objectEV{sub_i, 1}(end, :);
        sub_vv_coord = [structure_model_EV(1, 1:12)', structure_model_EV(1, 13:24)'];

        ctx_dist = [zeros(12, 1); ones(12, 1)];
        vv_grid = [[sub_vv_coord; sub_vv_coord], ctx_dist];
        vv_orth = [[[sub_vv_coord(:, 1), zeros(12, 1)]; [zeros(12, 1), sub_vv_coord(:, 2)]], ctx_dist];

        vv_para = [[[sub_vv_coord(:, 1), zeros(12, 1)]; [sub_vv_coord(:, 2), zeros(12, 1)]], ctx_dist];
        vv_EV01 = [[sub_vv_coord(:, 1);sub_vv_coord(:, 1)], ctx_dist];
        vv_EV02 = [[sub_vv_coord(:, 2);sub_vv_coord(:, 2)], ctx_dist];

        vv_rdm(:, 1) = zscore(pdist(vv_grid))';
        vv_rdm(:, 2) = zscore(pdist(vv_orth))';
        vv_rdm(:, 3) = zscore(pdist(vv_para))';

        roi_rdm = zscore(roi_cc_rdm{sub_i, roi_num, 4}');
        [b,bint,r] = regress(roi_rdm, [vv_rdm, ones(276, 1)]);
        sub_roi_beta(:, sub_i) = b;
    end

    disp(['figure : fitting the EVs with different structures to the vmPFC'])
    disp(' ')
    figure;
    group_color = {[255, 128, 0]./256, [7, 126, 151]./256};
    group_name = {'VTC_depended_group', 'Random_distribution_group'};
    bn = 1;
    bar_name = {'Grid', 'Orthogonal', 'Parallel'};
    for bi = 1:3
        for gl = 1:2
            gm_beta = sub_roi_beta(bi, group_label==gl);
            gm_beta = gm_beta(~isnan(gm_beta));
            gm_means = mean(gm_beta);
            gm_sems  = std(gm_beta)/sqrt(size(gm_beta, 2));
            % jitter = 0.2;
            b = bar(bn-0.5, gm_means, 0.6, 'FaceColor','none', 'EdgeColor','none');
            hold on
            % xi = bn-0.5 + (rand(size(gm_beta)) - 0.5) * 2 * jitter;
            % scatter(xi, gm_beta, 50, ones(1, 3).*0.2, 'filled', ...
            %     'MarkerFaceAlpha',0.7, 'MarkerEdgeAlpha',0.2);
            b.EdgeColor = group_color{gl};
            b.LineWidth = 5;
            errorbar(bn-0.5, gm_means, gm_sems, 'color', group_color{gl}, 'LineWidth',2, 'CapSize',12, 'LineStyle','none');
            bn = bn+1;
        end
    end
    box off
    set(gca, 'FontSize', 15, 'FontWeight', 'bold')
    set(gca,'LineWidth', 2)
    ylim([-0.5 1])
    xticks(1:2:6);
    xticklabels(bar_name);
    for gl = 1:2
        disp([group_name{gl}, ':'])
        for bi = 1:3
            group_beta = sub_roi_beta(bi, group_label==gl)';
            [h,p,ci,stats] = ttest(group_beta);
            disp(['ttest, ', bar_name{bi}, '~=0: p-value=', num2str(p)])
            tvalue(gl, bi) = stats.tstat;
        end
        disp(' ')
    end


    % VTC features represented in vmPFC on day 01
    cc_angle = (15:30:360)';
    cc_coord = [cosd(cc_angle), sind(cc_angle)];
    VTC_cc_rdm(:, 1) = zscore(pdist(cc_coord(:, 1)))';
    VTC_cc_rdm(:, 2) = zscore(pdist(cc_coord(:, 2)))';

    roi_num = 2;
    day01_vc_beta = nan(27, 5);
    for sub_i = 1:27
        structure_model_EV(1, :) = sub_structure_model_objectEV{sub_i, 1}(end, :);
        sub_vv_coord = [structure_model_EV(1, 1:12)', structure_model_EV(1, 13:24)'];
        ctx_vv_rdm(:, 1) = zscore(pdist(sub_vv_coord(:, 1)))';
        ctx_vv_rdm(:, 2) = zscore(pdist(sub_vv_coord(:, 2)))';

        roi_rdm = zscore(roi_cc_rdm{sub_i, roi_num, 1}');
        if ~isempty(roi_rdm)
            day01_vc_beta(sub_i, :) = regress(roi_rdm, [ctx_vv_rdm, VTC_cc_rdm, ones(66, 1)]);
        end
    end

    disp(['figure : Compare the representations of EVs and VTC features in the vmPFC on Day 1'])
    disp(' ')
    figure;
    group_color = {[255, 128, 0]./256, [7, 126, 151]./256};
    group_name = {'VTC_depended_group', 'Random_distribution_group'};
    bn = 1;
    bar_name = {'EV01', 'EV02', 'VTC f01', 'VTC f02'};
    for bi = 1:4
        for gl = 1:2
            gm_beta = day01_vc_beta(group_label==gl, bi)';
            gm_beta = gm_beta(~isnan(gm_beta));
            gm_means = mean(gm_beta);
            gm_sems  = std(gm_beta)/sqrt(size(gm_beta, 2));
            % jitter = 0.2;
            b = bar(bn-0.5, gm_means, 0.6, 'FaceColor','none', 'EdgeColor','none');
            hold on
            % xi = bn-0.5 + (rand(size(gm_beta)) - 0.5) * 2 * jitter;
            % scatter(xi, gm_beta, 50, ones(1, 3).*0.2, 'filled', ...
            %     'MarkerFaceAlpha',0.7, 'MarkerEdgeAlpha',0.2);
            b.EdgeColor = group_color{gl};
            b.LineWidth = 5;
            errorbar(bn-0.5, gm_means, gm_sems, 'color', group_color{gl}, 'LineWidth',2, 'CapSize',12, 'LineStyle','none');
            bn = bn+1;
        end
    end
    box off
    set(gca, 'FontSize', 15, 'FontWeight', 'bold')
    set(gca,'LineWidth', 2)
    ylim([-0.2 0.3])
    xticks(1:2:8);
    xticklabels(bar_name);
    for gl = 1:2
        disp([group_name{gl}, ':'])
        for bi = 1:4
            [h, p,ci,stats] = ttest(day01_vc_beta(group_label==gl, bi), 0, 'tail', 'right');
            disp(['ttest, ', bar_name{bi}, '>0: p-value=', num2str(p)])
            tvalue(gl, bi) = stats.tstat;
        end
        disp(' ')
    end



    % VTC features represented in vmPFC on day 03 (Regress the EV out)
    roi_num = 2;
    day03_ctx_vc_beta = nan(27, 2, 3);
    for sub_i = 1:27
        structure_model_EV(1, :) = sub_structure_model_objectEV{sub_i, 1}(end, :);
        sub_vv_coord = [structure_model_EV(1, 1:12)', structure_model_EV(1, 13:24)'];
        ctx_vv_rdm(:, 1) = zscore(pdist(sub_vv_coord(:, 1)))';
        ctx_vv_rdm(:, 2) = zscore(pdist(sub_vv_coord(:, 2)))';

        for ctx_i = 2:3
            roi_rdm = zscore(roi_cc_rdm{sub_i, roi_num, ctx_i}');
            [b,bint,r] = regress(roi_rdm, [ctx_vv_rdm, ones(66, 1)]);
            day03_ctx_vc_beta(sub_i, ctx_i-1, :) = regress(r, [VTC_cc_rdm, ones(66, 1)]);
        end
    end

    disp(['figure : Representations of VTC features in the vmPFC on Day 3'])
    disp(' ')
    figure;
    group_color = {[255, 128, 0]./256, [7, 126, 151]./256};
    group_name = {'VTC_depended_group', 'Random_distribution_group'};
    bn = 1;
    bar_name = {'CtxA VTC f01', 'CtxA VTC f02', 'CtxB VTC f01', 'CtxB VTC f02'};
    bar_name = [bar_name, bar_name];
    for gl = 1:2
        for ctx_i = 1:2
            for bi = 1:2
                gm_beta = day03_ctx_vc_beta(group_label==gl, ctx_i, bi)';
                gm_beta = gm_beta(~isnan(gm_beta));
                gm_means = mean(gm_beta);
                gm_sems  = std(gm_beta)/sqrt(size(gm_beta, 2));
                % jitter = 0.2;
                b = bar(bn-0.5, gm_means, 0.6, 'FaceColor','none', 'EdgeColor','none');
                hold on
                % xi = bn-0.5 + (rand(size(gm_beta)) - 0.5) * 2 * jitter;
                % scatter(xi, gm_beta, 50, ones(1, 3).*0.2, 'filled', ...
                %     'MarkerFaceAlpha',0.7, 'MarkerEdgeAlpha',0.2);
                b.EdgeColor = group_color{gl};
                b.LineWidth = 5;
                errorbar(bn-0.5, gm_means, gm_sems, 'color', group_color{gl}, 'LineWidth',2, 'CapSize',12, 'LineStyle','none');
                bn = bn+1;
            end
        end
    end
    box off
    set(gca, 'FontSize', 15, 'FontWeight', 'bold')
    set(gca,'LineWidth', 2)
    ylim([-0.0 0.4])
    xticks(1:1:16);
    xticklabels(bar_name);

    bn = 1;
    for gl = 1:2
        disp([group_name{gl}, ':'])
        for ctx_i = 1:2
            for bi = 1:2
                group_vc_beta = squeeze(day03_ctx_vc_beta(group_label==gl, ctx_i, bi));
                [h,p,ci,stats] = ttest(group_vc_beta, 0, 'tail', 'right');
                disp(['ttest, ', bar_name{bn}, '>0: p-value=', num2str(p)])
                bn = bn + 1;
            end
        end
        disp(' ')
    end
  

    % VTC feature model comparison through RSA
    cc_angle = (15:30:360)';
    cc_coord = [cosd(cc_angle), sind(cc_angle)];
    cc_grid = [cc_coord; cc_coord];
    cc_orth = [[cc_coord(:, 1), zeros(12, 1)]; [zeros(12, 1), cc_coord(:, 2)]];
    ra = 90;
    rotation = [[cosd(ra); sind(ra)], [-sind(ra); cosd(ra)]];
    cc_para = [[cc_coord(:, 1), zeros(12, 1)]; [zeros(12, 1), cc_coord(:, 2)]*rotation];
    cc_rdm(:, 1) = zscore(pdist(cc_grid))';
    cc_rdm(:, 2) = zscore(pdist(cc_orth))';
    cc_rdm(:, 3) = zscore(pdist(cc_para))';
    roi_num = 2;
    sub_roi_beta = [];
    for sub_i = 1:27
        target_rdm = zscore(roi_cc_rdm{sub_i, roi_num, 4}');
        [b,bint,r] = regress(target_rdm, [cc_rdm, ones(276, 1)]);
        sub_roi_beta(:, sub_i) = b;
    end
    disp(['figure : fitting the VTC features with different structures to the vmPFC'])
    disp(' ')
    figure;
    group_color = {[255, 128, 0]./256, [7, 126, 151]./256};
    group_name = {'VTC_depended_group', 'Random_distribution_group'};
    bn = 1;
    bar_name = {'Grid', 'Orthogonal', 'Parallel'};
    for bi = 1:3
        for gl = 1:2
            gm_beta = sub_roi_beta(bi, group_label==gl);
            gm_beta = gm_beta(~isnan(gm_beta));
            gm_means = mean(gm_beta);
            gm_sems  = std(gm_beta)/sqrt(size(gm_beta, 2));
            % jitter = 0.2;
            b = bar(bn-0.5, gm_means, 0.6, 'FaceColor','none', 'EdgeColor','none');
            hold on
            % xi = bn-0.5 + (rand(size(gm_beta)) - 0.5) * 2 * jitter;
            % scatter(xi, gm_beta, 50, ones(1, 3).*0.2, 'filled', ...
            %     'MarkerFaceAlpha',0.7, 'MarkerEdgeAlpha',0.2);
            b.EdgeColor = group_color{gl};
            b.LineWidth = 5;
            errorbar(bn-0.5, gm_means, gm_sems, 'color', group_color{gl}, 'LineWidth',2, 'CapSize',12, 'LineStyle','none');
            bn = bn+1;
        end
    end
    box off
    set(gca, 'FontSize', 15, 'FontWeight', 'bold')
    set(gca,'LineWidth', 2)
    ylim([-0.2 0.6])
    xticks(1:2:6);
    xticklabels(bar_name);
    for gl = 1:2
        disp([group_name{gl}, ':'])
        for bi = 1:3
            group_beta = sub_roi_beta(bi, group_label==gl)';
            [h,p,ci,stats] = ttest(group_beta);
            disp(['ttest, ', bar_name{bi}, '~=0: p-value=', num2str(p)])
            tvalue(gl, bi) = stats.tstat;
        end
        disp(' ')
    end




    outputArg1 = 'Analysis completed without error';


end

