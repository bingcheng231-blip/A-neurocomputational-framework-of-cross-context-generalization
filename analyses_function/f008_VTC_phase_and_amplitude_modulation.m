function [outputArg1] = f008_VTC_phase_and_amplitude_modulation(Re_analysis_data_path,exp_label)
%F008_VTC_PHASE_AND_AMPLITUDE_MODULATION 此处显示有关此函数的摘要
%   此处显示详细说明
    [group_label, rotation_direction, ~, ~] = s002_grouplabel(Re_analysis_data_path, exp_label);
    meanfunc_ROIlabel = load([Re_analysis_data_path, 'd009_sub_VTC_label.mat']).meanfunc_ROIlabel;

    % loading the voxel-level FFT results
    load([Re_analysis_data_path, 'd008_fMRI_sub_VTCvoxel_phase_and_amplitude.mat']);
    % polarhistogram
    disp(' ')
    disp('Figure :　Phase modulation of VTC voxels ')
    group_name = {'VTC_depended_group', 'Random_distribution_group'};
    context_color{1} = [0.4, 0.4, 0.4]; context_color{2} = [1, 0.15, 0.85]; context_color{3} = [0.15, 0.85, 1];
    for gi = 1:2
        disp([group_name{gi}])
        vtc_theta = cell2mat(sub_VTCvoxel_phase(group_label==gi));
        figure
        for ctx_i = 1:3
            ctx_vtc_theta = vtc_theta(:, ctx_i);
            polarhistogram(ctx_vtc_theta,36, 'FaceColor', context_color{ctx_i}, 'FaceAlpha',.4, 'EdgeColor', context_color{ctx_i})
            hold on
        end
        set(gca, 'FontSize', 15, 'FontWeight', 'bold')
        set(gca,'LineWidth',2)
    end


    voxel_prefer_direction = nan(27, 3, 2);
    for sub_i = [1, 3:27]
        % VTC label
        for ctx_i = 1:3
            vtc_theta = sub_VTCvoxel_phase{sub_i}(:, ctx_i);
            hem_vtc_theta{1} = vtc_theta(1:sum(meanfunc_ROIlabel{sub_i, 1}, 'all'));
            hem_vtc_theta{2} = vtc_theta(1+sum(meanfunc_ROIlabel{sub_i, 1}, 'all'):end);
            for hem_i = 1:2
                [voxel_num,edges] = histcounts(hem_vtc_theta{hem_i},72);
                theta_e0 = [];
                for ei = 1:72
                    theta_e0(ei) = mean(edges(ei:ei+1));
                end
                [x_list] = s009_Gaussian_perfer_direction(voxel_num,  [] , theta_e0, 32);
                voxel_prefer_direction(sub_i, ctx_i, hem_i) = x_list(1, 3);
            end
        end
    end

    disp(' ')
    disp('Figure :　Overall preferred direction of VTC')
    disp([group_name{1}])
    group_color{1, 1} = [255, 96, 0] ./255;   group_color{1, 2} = [255, 192, 128]./255;
    group_color{2, 1} = [0, 114, 189] ./255;  group_color{2, 2} = [64, 178, 254]./255;
    figure
    group_pd = squeeze(voxel_prefer_direction(group_label==1, :, 1));
    group_pd = rad2deg(group_pd);
    group_pd(group_pd<0) = abs(group_pd(group_pd<0));
    group_pd(group_pd>90) = 180-group_pd(group_pd>90);
    hem_pd{1, 1} = group_pd;
    vs = violinplot(group_pd,{'day01', 'day03_ctxA','day03_ctxB'},'ViolinColor',group_color{1, 1},'ShowData',true,'ViolinAlpha',0.5,'HalfViolin','left','BoxColor',[0 0 0],'MedianColor',[0.75 0.75 0.75]);
    hold on
    group_pd = squeeze(voxel_prefer_direction(group_label==1, :, 2));
    group_pd = rad2deg(group_pd);
    group_pd(group_pd<0) = abs(group_pd(group_pd<0));
    group_pd(group_pd>90) = 180-group_pd(group_pd>90);
    hem_pd{2, 1} = group_pd;
    vs = violinplot(group_pd,{'day01', 'day03 ctxA','day03 ctxB'},'ViolinColor',group_color{1, 2},'ShowData',true,'ViolinAlpha',0.5,'HalfViolin','right','BoxColor',[0.5 0.5 0.5],'MedianColor',[1 1 1]);
    box off
    set(gca, 'FontSize', 15, 'FontWeight', 'bold')
    set(gca,'LineWidth',2)
 
    VTC_pd = cell2mat(hem_pd);
    [p,h,stats] = signrank(VTC_pd(:, 2), VTC_pd(:, 3));
    disp(['day03 ctxA ~= day03 ctxB, signrank: p-value=', num2str(p)])


    figure
    disp([group_name{2}])
    group_pd = squeeze(voxel_prefer_direction(group_label==2, :, 1));
    group_pd = rad2deg(group_pd);
    group_pd(group_pd<0) = abs(group_pd(group_pd<0));
    group_pd(group_pd>90) = 180-group_pd(group_pd>90);
    hem_pd{1, 1} = group_pd;
    vs = violinplot(group_pd,{'day01', 'day03 ctxA','day03 ctxB'},'ViolinColor',group_color{2, 1},'ShowData',true,'ViolinAlpha',0.5,'HalfViolin','left','BoxColor',[0.5 0.5 0.5],'MedianColor',[1 1 1]);
    hold on
    group_pd = squeeze(voxel_prefer_direction(group_label==2, :, 2));
    group_pd = rad2deg(group_pd);
    group_pd(group_pd<0) = abs(group_pd(group_pd<0));
    group_pd(group_pd>90) = 180-group_pd(group_pd>90);
    hem_pd{2, 1} = group_pd;
    vs = violinplot(group_pd,{'day01', 'day03 ctxA','day03 ctxB'},'ViolinColor',group_color{2, 2},'ShowData',true,'ViolinAlpha',0.5,'HalfViolin','right','BoxColor',[0.5 0.5 0.5],'MedianColor',[1 1 1]);
    box off
    set(gca, 'FontSize', 15, 'FontWeight', 'bold')
    set(gca,'LineWidth',2)

    VTC_pd = cell2mat(hem_pd);
    [p,h,stats] = signrank(VTC_pd(:, 2), VTC_pd(:, 3));
    disp(['day03 ctxA ~= day03 ctxB, signrank: p-value=', num2str(p)])


    % the main theta change between context
    % day03 - day01
    disp(' ')
    disp('Figure :　Voxel phase: Day01 - Day03')
    context_color{1} = [0.4, 0.4, 0.4]; context_color{2} = [1, 0.15, 0.85]; context_color{3} = [0.15, 0.85, 1];
    for gi = 1:2
        disp([group_name{gi}])
        vtc_theta = cell2mat(sub_VTCvoxel_phase(group_label==gi));
        figure
        for ctx_i = 2:3
            ctx_vtc_theta = vtc_theta(:, ctx_i)-vtc_theta(:, 1);
            ctx_vtc_theta(ctx_vtc_theta>pi) = ctx_vtc_theta(ctx_vtc_theta>pi)-2*pi;
            ctx_vtc_theta(ctx_vtc_theta<-pi) = ctx_vtc_theta(ctx_vtc_theta<-pi)+2*pi;
            polarhistogram(ctx_vtc_theta,36, 'FaceColor', context_color{ctx_i}, 'FaceAlpha',.4, 'EdgeColor', context_color{ctx_i})
            hold on
        end
        box off
        set(gca, 'FontSize', 15, 'FontWeight', 'bold')
        set(gca,'LineWidth',2)
    end

    % day03 ctxA - ctxB
    disp(' ')
    disp('Figure : Voxel phase difference between　day03 context A and B')
    group_color = {[255, 128, 0]./256, [7, 126, 151]./256};
    figure
    for gi = 2:-1:1
        vtc_theta = cell2mat(sub_VTCvoxel_phase(group_label==gi));
        ctx_vtc_theta = vtc_theta(:, 2)-vtc_theta(:, 3);
        ctx_vtc_theta(ctx_vtc_theta>pi) = ctx_vtc_theta(ctx_vtc_theta>pi)-2*pi;
        ctx_vtc_theta(ctx_vtc_theta<-pi) = ctx_vtc_theta(ctx_vtc_theta<-pi)+2*pi;
        polarhistogram(ctx_vtc_theta,36, 'FaceColor', group_color{gi}, 'FaceAlpha',.4, 'EdgeColor', group_color{gi})
        hold on
    end
    box off
    set(gca, 'FontSize', 15, 'FontWeight', 'bold')
    set(gca,'LineWidth',2)



    % two VTC subgroup  day03 - day01
    disp(' ')
    disp('figure : two VTC subgroup  day03 - day01')
    rd_name = {'VTC_depended_clockwise', 'VTC_depended_anticlockwise'};
    context_color{1} = [0.4, 0.4, 0.4]; context_color{2} = [1, 0.15, 0.85]; context_color{3} = [0.15, 0.85, 1];
    rotation_list = [-1, 1];
    for ri = 1:2
        vtc_theta = cell2mat(sub_VTCvoxel_phase(rotation_direction==rotation_list(ri)));
        figure
        disp(rd_name{ri})
        for ctx_i = 2:3
            ctx_vtc_theta = vtc_theta(:, ctx_i)-vtc_theta(:, 1);
            ctx_vtc_theta(ctx_vtc_theta>pi) = ctx_vtc_theta(ctx_vtc_theta>pi)-2*pi;
            ctx_vtc_theta(ctx_vtc_theta<-pi) = ctx_vtc_theta(ctx_vtc_theta<-pi)+2*pi;
            polarhistogram(ctx_vtc_theta,36, 'FaceColor', context_color{ctx_i}, 'FaceAlpha',.4, 'EdgeColor', context_color{ctx_i})
            hold on
        end
        box off
        set(gca, 'FontSize', 15, 'FontWeight', 'bold')
        set(gca,'LineWidth',2)
    end


    % two VTC subgroup  day03 ctxA - day03 ctxB
    disp(' ')
    disp('figure : two VTC subgroup  day03 ctxA - day03 ctxB')
    rd_name = {'VTC_depended_clockwise', 'VTC_depended_anticlockwise'};
    rotation_list = [-1, 1];
    rotation_color = {[255, 216, 178]./256, [255, 164, 164]./256};
    figure
    for ri = 1:2
        vtc_theta = cell2mat(sub_VTCvoxel_phase(rotation_direction==rotation_list(ri)));
        ctx_vtc_theta = vtc_theta(:, 2)-vtc_theta(:, 3);
        ctx_vtc_theta(ctx_vtc_theta>pi) = ctx_vtc_theta(ctx_vtc_theta>pi)-2*pi;
        ctx_vtc_theta(ctx_vtc_theta<-pi) = ctx_vtc_theta(ctx_vtc_theta<-pi)+2*pi;
        polarhistogram(ctx_vtc_theta,36, 'FaceColor', rotation_color{ri}, 'FaceAlpha',.4, 'EdgeColor', rotation_color{ri})

        % histogram(ctx_vtc_theta,36, 'FaceColor', rotation_color{ri}, 'FaceAlpha',.4, 'EdgeColor', rotation_color{ri})
        hold on
    end
    box off
    set(gca, 'FontSize', 15, 'FontWeight', 'bold')
    set(gca,'LineWidth',2)



    % amplitude modulation at VTC voxel level
    for sub_i = [1, 3:27]
        VTCvoxel_phase = sub_VTCvoxel_phase{sub_i};
        VTCvoxel_cc_resp =  sub_VTCvoxel_amplitude{sub_i};
        for ctx_i = 1:3
            for bin_i = 1:36
                lb_theta = -pi + deg2rad((bin_i-1).*10);
                ub_theta = -pi + deg2rad(bin_i.*10);

                theta_bin = VTCvoxel_phase(:, ctx_i)>lb_theta & VTCvoxel_phase(:, ctx_i)<ub_theta;
                bin_amplitude(bin_i, ctx_i, sub_i) = mean(abs(VTCvoxel_cc_resp(theta_bin, ctx_i)));
            end
        end
    end
    disp(' ')
    disp('figure : amplitude modulation')
    context_color{1} = [0.4, 0.4, 0.4]; context_color{2} = [1, 0.15, 0.85]; context_color{3} = [0.15, 0.85, 1];
    for gl = 1:2
        figure
        disp(group_name{gl})
        theta = deg2rad(1:361);
        group_color = {[255, 128, 0]./256, [7, 126, 151]./256};
        ShadowAlpha = 0.2;
        background_MaxCorr = ceil(max(mean(squeeze(bin_amplitude(:, 2, :)), 2))*1.2)+1;
        background_LineWidth = 1.2;
        background_Linecolor = 0.7.*ones(1, 3);
        for ri = background_MaxCorr/4:background_MaxCorr/4:background_MaxCorr
            [x,y] = pol2cart(theta, ri);
            plot(x, y, 'LineWidth', background_LineWidth, 'color',  background_Linecolor)
            hold on
        end
        for ti = 0:30:330
            [x,y] = pol2cart(deg2rad(ti), (0:0.01:background_MaxCorr)');
            plot(x, y, 'LineWidth', background_LineWidth, 'color',  background_Linecolor)
            hold on
        end
        theta = deg2rad(-175:10:185);
        for ctx_i = 1:3
            group_bin_am = squeeze(bin_amplitude(:, ctx_i, group_label==gl))';
            group_bin_am(:, 37) = group_bin_am(:, 1);
            SEM=std(group_bin_am)./sqrt(sum(group_label==gl));
            [x1_SEM,y1_SEM] = pol2cart(theta, mean(group_bin_am)-SEM);
            [x2_SEM,y2_SEM] = pol2cart(theta, mean(group_bin_am)+SEM);
            fill([x1_SEM, fliplr(x2_SEM)], [y1_SEM, fliplr(y2_SEM)], 'k', 'FaceAlpha', ShadowAlpha, 'EdgeColor', 'none');
            hold on
            [x,y] = pol2cart(theta, mean(group_bin_am));
            plot(x, y, 'LineWidth', 1.7, 'color', context_color{ctx_i})
            hold on
        end
        axis equal
        axis off
     end


     figure
     disp(' ')
     disp('figure : 0°/180° direction voxel activity')
     bn = 1;
     bin_name = {'day01', 'day03 ctxA', 'day03 ctxB'};
     bin_name = [bin_name, bin_name];
     for gl = 1:2
         bin_num = [1, 36, 18, 19];
         group_bin_mean = squeeze(mean(bin_amplitude(bin_num, :, group_label==gl), 1));
         [p,h,stats] = signrank(group_bin_mean(2, :)', group_bin_mean(3, :)', 'tail', 'right');
         disp([group_name{gl}, ': day03 ctxA > day03 ctxB, signrank, p-value=', num2str(p)])
         for mi = 1:size(group_bin_mean, 1)
             gm_beta = group_bin_mean(mi, :);
             gm_means = mean(gm_beta);
             gm_sems  = std(gm_beta)/sqrt(length(gm_beta));
             % jitter = 0.2;
             b = bar(bn-0.5, gm_means, 0.6, 'FaceColor','none', 'EdgeColor','none');
             hold on
             % xi = mi-0.5 + (rand(size(gm_NeglogLik)) - 0.5) * 2 * jitter;
             % scatter(xi, gm_NeglogLik, 50, ones(1, 3).*0.2, 'filled', ...
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
     xticks(0.5:7);
     xticklabels(bin_name);

     figure
     disp(' ')
     disp('figure : 90°/270° direction voxel activity')
     bn = 1;
     bin_name = {'day01', 'day03 ctxA', 'day03 ctxB'};
     bin_name = [bin_name, bin_name];
     for gl = 1:2
         bin_num = [9, 10, 27, 28];
         group_bin_mean = squeeze(mean(bin_amplitude(bin_num, :, group_label==gl), 1));
         [p,h,stats] = signrank(group_bin_mean(2, :)', group_bin_mean(3, :)', 'tail', 'left');
         disp([group_name{gl}, ': day03 ctxA < day03 ctxB, signrank, p-value=', num2str(p)])
         for mi = 1:size(group_bin_mean, 1)
             gm_beta = group_bin_mean(mi, :);
             gm_means = mean(gm_beta);
             gm_sems  = std(gm_beta)/sqrt(length(gm_beta));
             % jitter = 0.2;
             b = bar(bn-0.5, gm_means, 0.6, 'FaceColor','none', 'EdgeColor','none');
             hold on
             % xi = mi-0.5 + (rand(size(gm_NeglogLik)) - 0.5) * 2 * jitter;
             % scatter(xi, gm_NeglogLik, 50, ones(1, 3).*0.2, 'filled', ...
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
     xticks(0.5:7);
     xticklabels(bin_name);



     % the voxel amplitude change between context
     % day03 - day01
     disp(' ')
     disp('figure : Voxel amplitude: Day01 - Day03 context')
     context_color{1} = [0.4, 0.4, 0.4]; context_color{2} = [1, 0.15, 0.85]; context_color{3} = [0.15, 0.85, 1];
     for gi = 1:2
         VTCvoxel_cc_resp =  cell2mat(sub_VTCvoxel_amplitude(group_label==gi));
         figure
         disp(group_name{gi})
         for ctx_i = 2:3
             ctx_vtc_resp = VTCvoxel_cc_resp(:, ctx_i)-VTCvoxel_cc_resp(:, 1);
             histogram(ctx_vtc_resp,36, 'FaceColor', context_color{ctx_i}, 'FaceAlpha',.4, 'EdgeColor', context_color{ctx_i})
             hold on
         end
         box off
         if gi == 1
             xlim([-30 60])
         else
             xlim([-30 50])
         end
         set(gca, 'FontSize', 15, 'FontWeight', 'bold')
         set(gca,'LineWidth',2)
     end


     % day03 ctxA - ctxB
     group_color = {[255, 128, 0]./256, [7, 126, 151]./256};
     figure
     disp(' ')
     disp('figure : Voxel amplitude: day03 ctxA - ctxB')
     for gi = 2:-1:1
         disp(group_name{gi})
         VTCvoxel_cc_resp =  cell2mat(sub_VTCvoxel_amplitude(group_label==gi));
         ctx_vtc_resp = VTCvoxel_cc_resp(:, 2)-VTCvoxel_cc_resp(:, 3);
         histogram(ctx_vtc_resp,36, 'FaceColor', group_color{gi}, 'FaceAlpha',.4, 'EdgeColor', group_color{gi})
         hold on
     end
     box off
     set(gca, 'FontSize', 15, 'FontWeight', 'bold')
     set(gca,'LineWidth',2)



     % two VTC subgroup  day03 ctxA - day03 ctxB
     disp(' ')
     disp('figure : two VTC subgroup  day03 ctxA - day03 ctxB')
     rd_name = {'VTC_depended_clockwise', 'VTC_depended_anticlockwise'};
     rotation_list = [-1, 1];
     rotation_color = {[255, 216, 178]./256, [255, 164, 164]./256};
     figure
     for ri = 1:2
         VTCvoxel_cc_resp = cell2mat(sub_VTCvoxel_amplitude(rotation_direction==rotation_list(ri)));
         ctx_vtc_resp = VTCvoxel_cc_resp(:, 2)-VTCvoxel_cc_resp(:, 3);
         histogram(ctx_vtc_resp,36, 'FaceColor', rotation_color{ri}, 'FaceAlpha',.4, 'EdgeColor', rotation_color{ri})
         hold on
     end
     box off
     set(gca, 'FontSize', 15, 'FontWeight', 'bold')
     set(gca,'LineWidth',2)

     outputArg1 = 'Analysis completed without error';
end


