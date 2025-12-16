function [outputArg1] = f011_Gaussian_MLP_simulation(Re_analysis_data_path,exp_label, re_analysis)
%F011_GENERATE_MLP_SIMULATION 此处显示有关此函数的摘要
%   此处显示详细说明
    [group_label, rotation_direction, sub_reward, ~] = s002_grouplabel(Re_analysis_data_path, exp_label);
    sub_mlp = load([Re_analysis_data_path, 'd011_mGaussianMLP_WI0.01_phaseR3_amplitudeR1.mat']).sub_mlp;
    it4=wd_compit4_1;

    % training loss
    disp(' ')
    disp('figure : training loss')
    delta_t = 100;
    sub_el = [];
    for mi  = 1:4
        for sub_i = 1:size(sub_mlp, 1)
            sub_echo_loss = sub_mlp{sub_i, 6, mi};
            for ei = 1:size(sub_mlp{sub_i, 6, mi}(1, :), 2)./delta_t
                % mean loss
                sub_el(sub_i, ei, mi) = mean(sub_echo_loss(1+delta_t*(ei-1):delta_t+delta_t*(ei-1)));
            end
        end
    end
    group_color = {[255, 128, 0]./256, [7, 126, 151]./256};
    x = linspace(1,size(sub_mlp{sub_i, 6, 1}(1, :), 2)*24,40);
    line_type = {'--', '-o', '-x', '-^'};
    % line type => ['none model', 'phase modulation model', 'amplitude modulation model', 'both modulation model']
    group_name = {'VTC_dependent', 'Random_distribution'};
    for gl = 1:2
        figure
        disp(group_name{gl})
        for mi = 1:4
            y = squeeze(sub_el(group_label==gl, :, mi));
            SEM=std(y)./sqrt(size(y, 1));
            s = shadedErrorBar(x,mean(y, 1),SEM, 'lineProps', '-');
            s.patch.FaceColor = [0.5, 0.5, 0.5];
            set(s.edge,'LineWidth',0.5, 'Color', [0.7, 0.7, 0.7])
            hold on
            plot(x, mean(y, 1), line_type{mi}, 'LineWidth', 1.5, 'color',  group_color{gl})
            hold on
        end
        ylim([0, 60])
        set(gca, 'FontSize', 15, 'FontWeight', 'bold')
        set(gca,'LineWidth',2)
    end



    % Firing rate map of amplitude-modulated VTCneurons
    disp(' ')
    disp('figure : Firing rate map of amplitude-modulated VTCneurons')
    [X, Y] = meshgrid(-1.5:0.1:1.5);
    Y=-Y;
    cc_angle = (15:30:360)';
    VTCneuron_mu = [cosd(cc_angle), sind(cc_angle)];
    VTC_fr_max = 0.7879;
    pos_num = numel(X); neuron_num = 12;
    sub_activation_map = [];
    for sub_i = 1:36
        mlp_amplitude = squeeze(sub_mlp{sub_i, 7, 3}(40, :, :));
        for ctx_i = 1:2
            vtc_sigm = zeros(2, 2);
            activation_map = nan(pos_num, neuron_num);
            for vtc_ni = 1:12
                vtc_sigma(1, 1) = 0.202 + abs(VTCneuron_mu(vtc_ni, 1).*mlp_amplitude(vtc_ni, ctx_i));
                vtc_sigma(2, 2) = 0.202 + abs(VTCneuron_mu(vtc_ni, 2).*mlp_amplitude(vtc_ni, ctx_i));
    
                fr = mvnpdf([X(:),Y(:)],VTCneuron_mu(vtc_ni, :),vtc_sigma)./VTC_fr_max;
                activation_map(:, vtc_ni) = fr;
            end
            sub_activation_map(:, :, ctx_i, sub_i) = activation_map;
        end
    end
    group_name = {'VTC_dependent', 'Random_distribution'};
    ctx_name = {'contextA', 'contextB'};
    for ctx_i = 1:2
        for gl = 1:2
            group_mean_activation_map = mean(squeeze(sub_activation_map(:, :, ctx_i, group_label==gl)), 3);
            disp([ctx_name{ctx_i}, ' ', group_name{gl}])
            figure
            for ni = 1:neuron_num
                fr = reshape(group_mean_activation_map(:, ni),size(X));
                aps = fr;
                map_corr=[];
                for i=1:size(aps,1)
                    for j=1:size(aps,2)
                        map_corr(i,j,1:3)=it4(max(min(round(255*aps(i,j))+1,256),1),1:3);
                    end
                end
                subplot(3, 4, ni)
                imshow(map_corr)
                shading interp
            end
        end
    end

  
    % Shifting of phase-modulated VTCneurons
    % two groups
    disp(' ')
    disp('figure : Shifting of phase-modulated VTCneurons')
    ctx_mark = {"o", '^'};
    group_name = {'VTC_dependent', 'Random_distribution'};
    center_colormap = hsv(360);
    center_color = center_colormap(cc_angle, :);
    for gl = 1:2
        disp(group_name{gl})
        figure
        for ctx_i = 1:2
            phase = [];
            for sub_i = 1:36
                phase(sub_i, :) = deg2rad(cc_angle)+sub_mlp{sub_i, 8, 2}(40, :, ctx_i)';
            end
            group_mean_phase = mean(phase(group_label==gl, :))';
            VTCneuron_mu = [cos(group_mean_phase), sin(group_mean_phase)];  
            % show the firing rate map of VTC neurons
            scatter(VTCneuron_mu(:, 1), VTCneuron_mu(:, 2), 240, center_color, 'filled', ctx_mark{ctx_i})
            hold on
        end
        xlim([-1 1])
        ylim([-1 1])
        axis square
        axis off
    end
 
    % two subgroups
    disp(' ')
    ctx_mark = {"o", '^'};
    rd_list = [-1, 1];
    rd_name = {'VTC_depended_clockwise', 'VTC_depended_anticlockwise'};
    for gl = 1:2
        disp(rd_name{gl})
        figure
        for ctx_i = 1:2
            phase = [];
            for sub_i = 1:36
                phase(sub_i, :) = deg2rad(cc_angle)+sub_mlp{sub_i, 8, 2}(40, :, ctx_i)';
            end
            group_mean_phase = mean(phase(rotation_direction==rd_list(gl), :))';
            VTCneuron_mu = [cos(group_mean_phase), sin(group_mean_phase)];
            % show the firing rate map of VTC neurons
            scatter(VTCneuron_mu(:, 1), VTCneuron_mu(:, 2), 240, center_color, 'filled', ctx_mark{ctx_i})
            hold on
        end
        xlim([-1 1])
        ylim([-1 1])
        axis square
        axis off
    end




    %  Geometric relationships of VTC features across different contexts in MLP models.
    disp(' ')
    if re_analysis
        disp('Analysis geometric relationship of VTC features in MLP.... ')
        disp('This part of the analysis needs to run for a long time...')
        cc_angle = (15:30:360)';
        cc_coord = [cosd(cc_angle), sind(cc_angle)];
        rep_num = 32;
        LB = [0, zeros(1, 4)]; % Lower bound
        UB = [10, ones(1, 4)]; % Upper bound
        params_num = length(LB);
        OPTIM_options = optimset('Display', 'off') ;
        roi_vtc_features_dist = []; roi_cc_rotation_rr = [];
        for mi = 1:4
            for sub_i = 1:size(sub_dir, 1)
                vtc_representation =  sub_mlp{sub_i, 1, mi};
                for ei = 1:40
                    ctx_rdm = pdist(squeeze(vtc_representation(ei, :, :)))';
                    target_rdm = zscore(ctx_rdm);
                    angle_range = 1:3:361;
                    parfor ra_ii = 1:121
                        ra_i = angle_range(ra_ii);
                        rotation = [[cosd(ra_i); sind(ra_i)], [-sind(ra_i); cosd(ra_i)]];
                        X0 = LB + (UB-LB).*rand(rep_num, params_num); % random initialisation
                        theta_ML_all = nan(rep_num, params_num); NeglogLik_ML_all = nan(rep_num, 1);
                        for ri = 1:rep_num
                            lastwarn('');
                            [theta_ML_all(ri,:), NeglogLik_ML_all(ri)] = fmincon(@(theta) sumsqr(target_rdm'-s008_ROIs_rdm_ctxDist_calc(cc_coord, rotation, theta)),X0(ri,:),[],[],[],[],LB', UB',[],OPTIM_options);
                            [warnMsg, warnId] = lastwarn;
                            rep_c = 0;
                            while ~isempty(warnMsg)
                                disp('Recalculating...')
                                lastwarn('');
                                X0 = LB + (UB-LB).*rand(rep_num, params_num); % random initialisation
                                [theta_ML_all(ri,:), NeglogLik_ML_all(ri)] = fmincon(@(theta) sumsqr(target_rdm'-s008_ROIs_rdm_ctxDist_calc(cc_coord, rotation, theta)),X0(ri,:),[],[],[],[],LB', UB',[],OPTIM_options);
                                [warnMsg, warnId] = lastwarn;
                                rep_c = rep_c+1;
                                if rep_c>100
                                    error('Error initialization')
                                end
                            end
                        end
                        theta_all = [theta_ML_all, NeglogLik_ML_all];
                        theta_ML = sortrows(theta_all, size(theta_all, 2));
                        roi_vtc_features_dist(sub_i, mi, ei, ra_ii, :) = theta_ML(1, 1:end-1);
                        roi_cc_rotation_rr(sub_i, mi, ei, ra_ii) = corr(target_rdm, s008_ROIs_rdm_ctxDist_calc(cc_coord, rotation, theta_ML(1, 1:end-1))');
                    end
                end
            end
        end
    else
        roi_cc_rotation_rr = load([Re_analysis_data_path, 'd012_mGaussianMLP_WI0.01_phaseR3_amplitudeR1_VTCfeature_geomtery.mat']).roi_cc_rotation_rr;
        roi_vtc_features_dist = load([Re_analysis_data_path, 'd012_mGaussianMLP_WI0.01_phaseR3_amplitudeR1_VTCfeature_geomtery.mat']).roi_vtc_features_dist;
    end


    disp(' ')
    if re_analysis
        sub_model_ra = [];
        disp(' ')
        disp('Estimate the optimal feature rotation angle...')
        sa_theta = deg2rad(-179:3:178)';
        for sub_i = 1:36
            tic
            for mi = 1:4
                for ei = 1:40
                    ra = squeeze(roi_cc_rotation_rr(sub_i, mi, ei, :));
                    sub_ra = [ra(61:120); ra(1:60)];
                    [x_list] = s009_Gaussian_perfer_direction(sub_ra,  [] , sa_theta, 32);
                    sub_model_ra(sub_i, mi, ei) = rad2deg(x_list(1, 3));
                end
            end
            t = toc;
            disp(['Sub_', num2str(sub_i), '_done_t', num2str(t)])
        end
    else
        sub_model_ra = load([Re_analysis_data_path, 'd013_MLP_theoretical_VTCfeatures_representation_VTClayer.mat']).sub_model_ra;
    end

    disp(' ')
    disp('figure : VTCfeatures rotation in the VTClayer of MLP during the training  ')
    % VTCfeature rotation
    rotation_angle = [];
    for sub_i = 1:36
        for mi = 1:4
            for ei = 1:40
                rotation_angle(sub_i, ei, mi) = sub_model_ra(sub_i, mi, ei);
            end
        end
    end
    show_angle = abs(rotation_angle);
    group_color = {[255, 128, 0]./256, [7, 126, 151]./256};
    x = linspace(1,4000*24,40);
    line_type = {'--', '-o', '-x', '-^'}; % ['none model', 'phase modulation model', 'amplitude modulation model', 'both modulation model']
    group_name = {'VTC_dependent', 'Random_distribution'};
    for gl = 1:2
        figure
        disp(group_name{gl})
        for mi = 1:4
            y = squeeze(show_angle(group_label==gl, :, mi));
            SEM=std(y)./sqrt(size(y, 1));
            s = shadedErrorBar(x,mean(y, 1),SEM, 'lineProps', '-');
            s.patch.FaceColor = [0.5, 0.5, 0.5];
            set(s.edge,'LineWidth',0.5, 'Color', [0.7, 0.7, 0.7])
            hold on
            plot(x, mean(y, 1), line_type{mi}, 'LineWidth', 1.5, 'color',  group_color{gl})
            hold on
        end
        ylim([0, 80])
        set(gca, 'FontSize', 15, 'FontWeight', 'bold')
        set(gca,'LineWidth',2)
    end



    disp(' ')
    if re_analysis
        % theoretical VTCfeature representation
        disp(' ')
        disp('theoretical VTCfeature representation of Gaussian MLP')
        cc_angle = (15:30:360)';
        cc_coord = [cosd(cc_angle), sind(cc_angle)];
        rep_num = 32;
        LB = [0, zeros(1, 4)]; % Lower bound
        UB = [10, ones(1, 4)]; % Upper bound
        params_num = length(LB);
        OPTIM_options = optimset('Display', 'off') ;
        sub_tVTCfeatures = cell(36, 7, 4);
        for sub_i = 1:36
            tic
            for mi = 1:4
                vtc_layer_representation =  sub_mlp{sub_i, 1, mi};
                for ei = 1:40
                    ctx_rdm = pdist(squeeze(vtc_layer_representation(ei, :, :)))';
                    target_rdm = zscore(ctx_rdm);
                    ra_i = sub_model_ra(sub_i, mi, ei);
                    rotation = [[cosd(ra_i); sind(ra_i)], [-sind(ra_i); cosd(ra_i)]];
                    X0 = LB + (UB-LB).*rand(rep_num, params_num); % random initialisation
                    theta_ML_all = nan(rep_num, params_num); NeglogLik_ML_all = nan(rep_num, 1);
                    parfor ri = 1:rep_num
                        lastwarn('');
                        X0 = LB + (UB-LB).*rand(rep_num, params_num); % random initialisation
                        [theta_ML_all(ri,:), NeglogLik_ML_all(ri)] = fmincon(@(theta) sumsqr(target_rdm'-s008_ROIs_rdm_ctxDist_calc(cc_coord, rotation, theta)),X0(ri,:),[],[],[],[],LB', UB',[],OPTIM_options);
                        [warnMsg, warnId] = lastwarn;
                        rep_c = 0;
                        while ~isempty(warnMsg)
                            disp('Recalculating...')
                            lastwarn('');
                            X0 = LB + (UB-LB).*rand(rep_num, params_num); % random initialisation
                            [theta_ML_all(ri,:), NeglogLik_ML_all(ri)] = fmincon(@(theta) sumsqr(target_rdm'-s008_ROIs_rdm_ctxDist_calc(cc_coord, rotation, theta)),X0(ri,:),[],[],[],[],LB', UB',[],OPTIM_options);
                            [warnMsg, warnId] = lastwarn;
                            rep_c = rep_c+1;
                            if rep_c>100
                                error('Error initialization')
                            end
                        end
                    end

                    theta_all = [theta_ML_all, NeglogLik_ML_all];
                    theta_ML = sortrows(theta_all, size(theta_all, 2));
                    sub_tVTCfeatures{sub_i, mi, ei, 1} = theta_ML;
                    [model_rdm, ctx_vxy] = s008_ROIs_rdm_ctxDist_calc(cc_coord, rotation, theta_ML(1, 1:end-1));
                    sub_tVTCfeatures{sub_i, mi, ei, 2} = model_rdm';
                    sub_tVTCfeatures{sub_i, mi, ei, 3} = corr(target_rdm, sub_tVTCfeatures{sub_i, mi, ei, 2});
                    sub_tVTCfeatures{sub_i, mi, ei, 4} = ctx_vxy;
                end
            end
            t = toc;
            disp(['Sub_', num2str(sub_i, '%02d'), '_t', num2str(t)])
        end
    else
        sub_tVTCfeatures = load([Re_analysis_data_path, 'd013_MLP_theoretical_VTCfeatures_representation_VTClayer.mat']).sub_tVTCfeatures;
    end



    % task-irrelvent VTCfeatures compression
    disp(' ')
    disp('figure : task-irrelvent VTCfeatures compression')
    average_compression = [];
    for sub_i = 1:36
        for mi = 1:4
            for ei = 1:40
                theta_ML = sub_tVTCfeatures{sub_i, mi, ei, 1};
                cc_compression =  mean(theta_ML(1:10, 2:5));
                ctx_compression(1) = log(cc_compression(1)./cc_compression(2));
                ctx_compression(2) = log(cc_compression(4)./cc_compression(3));
                average_compression(sub_i, ei, mi) = mean(ctx_compression);
            end
        end
    end
    group_color = {[255, 128, 0]./256, [7, 126, 151]./256};
    x = linspace(1,4000*24,40);
    line_type = {'--', '-o', '-x', '-^'}; % ['none model', 'phase modulation model', 'amplitude modulation model', 'both modulation model']
    group_name = {'VTC_dependent', 'Random_distribution'};
    for gl = 1:2
        figure
        disp(group_name{gl})
        for mi = 1:4
            y = squeeze(average_compression(group_label==gl, :, mi));
            SEM=std(y)./sqrt(size(y, 1));
            s = shadedErrorBar(x,mean(y, 1),SEM, 'lineProps', '-');
            s.patch.FaceColor = [0.5, 0.5, 0.5];
            set(s.edge,'LineWidth',0.5, 'Color', [0.7, 0.7, 0.7])
            hold on
            plot(x, mean(y, 1), line_type{mi}, 'LineWidth', 1.5, 'color',  group_color{gl})
            hold on
        end
        ylim([-1, 0.5])

        set(gca, 'FontSize', 15, 'FontWeight', 'bold')
        set(gca,'LineWidth',2)
    end


    % VTC features' representational geometry in MLP
    disp(' ')
    disp('figure :  VTC features representational geometry in MLP ')
    model_name = {'none', 'phase', 'amplitude', 'both'};
    for model_num = 1:4
        disp(model_name{mi})
        model_cc_rotation_rr = squeeze(roi_cc_rotation_rr(:, :, 40, :));
        theta = deg2rad(1:3:361);
        figure
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
            rotation_rr = squeeze(model_cc_rotation_rr(rotation_direction==rd_list(rd_i), model_num, :));
            SEM=std(rotation_rr)./sqrt(sum(rotation_direction==rd_list(rd_i)));
            [x1_SEM,y1_SEM] = pol2cart(theta, mean(rotation_rr)-SEM);
            [x2_SEM,y2_SEM] = pol2cart(theta, mean(rotation_rr)+SEM);
            fill([x1_SEM, fliplr(x2_SEM)], [y1_SEM, fliplr(y2_SEM)], 'k', 'FaceAlpha', ShadowAlpha, 'EdgeColor', 'none');
            hold on
            [x,y] = pol2cart(theta, mean(rotation_rr));
            plot(x, y, 'LineWidth', 1.7, 'color', rotation_color{rd_i})
            hold on
        end
        rotation_rr = squeeze(model_cc_rotation_rr(group_label==2, model_num, :));
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
    end


    
    % individual level show case
    disp(' ')
    disp('figure : individual level show case of local phase modulation')
    disp('show category info: ')
    sub_i = 2;
    mi = 2;
    trial_reward = sub_reward{sub_i};
    catergoy_value = [trial_reward(:, 1); trial_reward(:, 2)]./10;
    catergoy_value = catergoy_value./2 + 0.5;
    aps = catergoy_value;
    value_corr=[];
    for i=1:size(aps,1)
        value_corr(i,1:3)=it4(max(min(round(255*aps(i))+1,256),1),1:3);
    end
    ctx_name = {'contextA', 'contextB'};
    for ctx_i = 1:2
        phase = [];
        phase(:, 1) = deg2rad(cc_angle)+sub_mlp{sub_i, 8, mi}(40, :, ctx_i)';
        VTCneuron_mu = [cos(phase), sin(phase)];
        % show the firing rate map of VTC neurons
        center_colormap = hsv(360);
        center_color = center_colormap(cc_angle, :);
        figure
        disp(ctx_name{ctx_i})
        scatter(VTCneuron_mu(:, 1), VTCneuron_mu(:, 2), 180, center_color, 'filled')
        hold on
        plot([0, 0], [-1.2, 1.2], '-k', 'LineWidth', 3)
        hold on
        plot([-1.2, 1.2], [0, 0], '-k', 'LineWidth', 3)
    
        xlim([-1.2 1.2])
        ylim([-1.2 1.2])
        axis square
        axis off
    end

    disp('show target rewards info: ')
    ctx_name = {'contextA', 'contextB'};
    for ctx_i = 1:2
        phase = [];
        phase(:, 1) = deg2rad(cc_angle)+sub_mlp{sub_i, 8, mi}(40, :, ctx_i)';
        VTCneuron_mu = [cos(phase), sin(phase)];
        % show the firing rate map of VTC neurons
        center_colormap = hsv(360);
        center_color = center_colormap(cc_angle, :);
        figure
        disp(ctx_name{ctx_i})
        scatter(VTCneuron_mu(:, 1), VTCneuron_mu(:, 2), 180, value_corr(1+(ctx_i-1)*12:12+(ctx_i-1)*12, :), 'filled')
        hold on
        plot([0, 0], [-1.2, 1.2], '-k', 'LineWidth', 3)
        hold on
        plot([-1.2, 1.2], [0, 0], '-k', 'LineWidth', 3)

        xlim([-1.2 1.2])
        ylim([-1.2 1.2])
        axis square
        axis off
    end

    outputArg1 = 'Analysis completed without error';
end

