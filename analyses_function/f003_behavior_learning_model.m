function [outputArg1] = f003_behavior_learning_model(Re_analysis_data_path, exp_label)
%F003_BEHAVIOR_LEARNING_MODEL 此处显示有关此函数的摘要
%   此处显示详细说明
    sub_learning_trials_info = load([Re_analysis_data_path, 'd004_Sub_learning_trials_info.mat']).sub_learning_trials_info;
    [group_label, ~, sub_traget_reward, exp_index] = s002_grouplabel(Re_analysis_data_path, exp_label);
    sub_learning_trials = sub_learning_trials_info(exp_index, :);
    it4=wd_compit4_1;
    close all

    all_behavior_model = cell(size(sub_learning_trials, 1), 2);

    % naive model
    disp(' Run naive model... ')
    rep_num = 256;
    lr_LB = 0.01;
    STRUCT_flag = false;
    decayed_lr_flag = true;
    prior_structure = [];
    for sub_i = 1:size(sub_learning_trials, 1)
        sub_choices = cell2mat(sub_learning_trials{sub_i, 2}');
        stimuli_id = (sub_choices(:, 2)-1)*12+sub_choices(:, 3);
        outcomes = sub_choices(:, 4);
        choices = sub_choices(:, 5);
       
        if decayed_lr_flag
            % [learning_rate, decay_rate, decay_steps, beta]
            LB = zeros(1, 4); % Lower bound 
            UB = [1, 1, 1080, 10.*ones(1, 1)]; % Upper bound 
            X0 = LB + (UB-LB).*rand(rep_num, size(LB, 2)); % random initialisation
        else
            % [learning_rate, beta]
            LB = zeros(1, 2); % Lower bound
            UB = [lr_LB, 10.*ones(1, 1)]; % Upper bound 
            X0 = LB + (UB-LB).*rand(rep_num, size(LB, 2)); % random initialisation
        end
        if STRUCT_flag
            LB = [LB, -1*ones(1, 276)]; % -1 is lower bound for cross-terms
            UB = [UB, ones(1, 276)]; % 1 is upper  bound for cross-terms
            X0 = [X0, -1 + 2*rand(rep_num,276)]; %  random initialisation
        end
    
        energy = @s003_ctx_RL_LogLikelihood;
        OPTIM_options = optimset('Display', 'off') ;
    
        NeglogLik_ML_all = nan(rep_num,1);
        theta_ML_all= nan(rep_num, size(LB, 2)); 
        parfor ri = 1:rep_num
            [theta_ML_all(ri,:), NeglogLik_ML_all(ri, 1)] = fmincon(@(theta) energy(theta,stimuli_id,outcomes,choices,STRUCT_flag,decayed_lr_flag, prior_structure),X0(ri,:)',[],[],[],[],LB', UB',[],OPTIM_options);
        end
        theta_all = [theta_ML_all, NeglogLik_ML_all];
        theta_ML = sortrows(theta_all, size(theta_all, 2));
        X0 = sortrows([X0, NeglogLik_ML_all], size(theta_all, 2));
        all_behavior_model{sub_i, 1} = theta_ML;
        t = toc;
        disp(['naive model: ', sub_learning_trials{sub_i, 1}, '_done_t', num2str(t)])
    end


    % structure model
    disp(' Run structure model... ')
    rep_num = 256;
    lr_LB = 0.01;
    STRUCT_flag = true;
    decayed_lr_flag = true;
    prior_structure = [];
    for sub_i = 1:size(sub_learning_trials, 1)
        sub_choices = cell2mat(sub_learning_trials{sub_i, 2}');
        stimuli_id = (sub_choices(:, 2)-1)*12+sub_choices(:, 3);
        outcomes = sub_choices(:, 4);
        choices = sub_choices(:, 5);
       
        if decayed_lr_flag
            % [learning_rate, decay_rate, decay_steps, beta]
            LB = zeros(1, 4); % Lower bound 
            UB = [1, 1, 1080, 10.*ones(1, 1)]; % Upper bound 
            X0 = LB + (UB-LB).*rand(rep_num, size(LB, 2)); % random initialisation
        else
            % [learning_rate, beta]
            LB = zeros(1, 2); % Lower bound
            UB = [lr_LB, 10.*ones(1, 1)]; % Upper bound 
            X0 = LB + (UB-LB).*rand(rep_num, size(LB, 2)); % random initialisation
        end
        if STRUCT_flag
            LB = [LB, -1*ones(1, 276)]; % -1 is lower bound for cross-terms
            UB = [UB, ones(1, 276)]; % 1 is upper  bound for cross-terms
            X0 = [X0, -1 + 2*rand(rep_num,276)]; %  random initialisation 
        end
        energy = @s003_ctx_RL_LogLikelihood;
        OPTIM_options = optimset('Display', 'off') ;
    
        NeglogLik_ML_all = nan(rep_num,1);
        theta_ML_all= nan(rep_num, size(LB, 2)); 
        parfor ri = 1:rep_num
            [theta_ML_all(ri,:), NeglogLik_ML_all(ri, 1)] = fmincon(@(theta) energy(theta,stimuli_id,outcomes,choices,STRUCT_flag,decayed_lr_flag, prior_structure),X0(ri,:)',[],[],[],[],LB', UB',[],OPTIM_options);
        end
        theta_all = [theta_ML_all, NeglogLik_ML_all];
        theta_ML = sortrows(theta_all, size(theta_all, 2));
        X0 = sortrows([X0, NeglogLik_ML_all], size(theta_all, 2));
        all_behavior_model{sub_i, 2} = theta_ML;
        t = toc;
        disp(['structure model: ', sub_learning_trials{sub_i, 1}, '_done_t', num2str(t)])
    end


    % behavior model comparison
    for sub_i = 1:size(sub_learning_trials, 1)
        model_theta = all_behavior_model{sub_i, 1};
        model_NeglogLik(sub_i, 1) = model_theta(1, end);

        model_theta = all_behavior_model{sub_i, 2};
        model_NeglogLik(sub_i, 2) = model_theta(1, end);
    end
    disp(' ')
    disp('behavior model comparison, signrank:  ')
    disp('NeglogLik: Navie > strcture')
    group_name = {'VTC-dependent', 'random distribution'};
    group_color = {[255, 128, 0]./256, [7, 126, 151]./256};
    disp('figure : model neglogLik of each group')
    for gl = 1:2
        group_NeglogLik = model_NeglogLik(group_label==gl, :);
        [p,h,stats] = signrank(group_NeglogLik(:, 1), group_NeglogLik(:, 2), 'tail', 'right');
        disp([group_name{gl}, ' p-value:', num2str(p)])

        figure;
        group_means = mean(group_NeglogLik);
        group_sems  = std(group_NeglogLik)/sqrt(size(group_NeglogLik, 1));
        jitter = 0.2;
        for bi = 1:2
            b = bar(bi-0.5, group_means(bi), 0.6, 'FaceColor','none', 'EdgeColor','none');
            hold on
            xi = bi-0.5 + (rand(size(group_NeglogLik(:, bi))) - 0.5) * 2 * jitter;
            scatter(xi, group_NeglogLik(:, bi), 50, ones(1, 3).*0.2, 'filled', ...
                'MarkerFaceAlpha',0.7, 'MarkerEdgeAlpha',0.2);
            b.EdgeColor = group_color{gl};   
            b.LineWidth = 5;
        end
        errorbar((1:2)-0.5, group_means, group_sems, 'color', group_color{gl}, 'LineWidth',3, 'CapSize',20, 'LineStyle','none');
        box off
        set(gca, 'FontSize', 15, 'FontWeight', 'bold')
        set(gca,'LineWidth', 2)
        ylim([0 800])
        xticks(0.5:3);
        xticklabels({'Naive','Structure'});
    end


    % Fitting different RDM models to the object relations derived from the structural model
    cc_angle = (15:30:360)';
    cc_coord = [cosd(cc_angle), sind(cc_angle)];
    cc_rdm = zscore(pdist([cc_coord;cc_coord]))'; % category similarity rdm
    mix_lb = nchoosek(1:3, 2);
    repeatability =[]; model_beta = [];
    for sub_i = 1:size(sub_learning_trials, 1)
        trial_reward = sub_traget_reward{sub_i, 1};
        % traget reward rdm
        reward_rdm = zscore(pdist([trial_reward(:, 1); trial_reward(:, 2)]))'; 
        % context info rdm
        ctx_rdm = zscore(pdist([zeros(12, 1); ones(12, 1)]))';
        % mix effect 
        model_variables = [cc_rdm, reward_rdm, ctx_rdm];
        mix_effect = nan(276, 4);
        for mi = 1:3
            mix_effect(:, mi) = zscore(model_variables(:, mix_lb(mi, 1)) .* model_variables(:, mix_lb(mi, 2)));
        end
        mix_effect(:, 4) = zscore(model_variables(:, 1) .* model_variables(:, 2) .* model_variables(:, 3));

        % object relations derived from the structural model
        model_theta = all_behavior_model{sub_i, 2};
        object_relations = model_theta(1:5, end-276:end-1)';
        object_relations = 1-object_relations;
        repeatability(sub_i, 1) = mean(1-pdist(object_relations', 'correlation'));

        % GLM
        model_effect = [model_variables, mix_effect, ones(276, 1)];
        % Using the mean object relations to improve robustness
        model_beta(sub_i, :) = lsqnonneg(model_effect, zscore(mean(object_relations, 2)));

        % show the mean object relations
        squared_object_relations = squareform(mean(object_relations, 2));
        estimation_structure(sub_i, :, :) = squared_object_relations;
        value_sequence = [sortrows([(1:12)', trial_reward(:, 1)], 2); sortrows([(13:24)', trial_reward(:, 2)], 2)];
        rearrangedMatrix = squared_object_relations(value_sequence(:, 1), value_sequence(:, 1));
        rearranged_structure(sub_i, :, :) = rearrangedMatrix;
    end


    % show the object relations repeatability
    disp(' ')
    disp('figure : the repeatability of the object relations derived from the structure model')
    figure;
    for gl = 1:2
        group_repeatability = repeatability(group_label==gl, 1);
        group_means = mean(group_repeatability);
        group_sems  = std(group_repeatability)/sqrt(size(group_repeatability, 1));
        jitter = 0.2;
        b = bar(gl-0.5, group_means, 0.6, 'FaceColor','none', 'EdgeColor','none');
        hold on
        xi = gl-0.5 + (rand(size(group_repeatability)) - 0.5) * 2 * jitter;
        scatter(xi, group_repeatability, 50, ones(1, 3).*0.2, 'filled', ...
            'MarkerFaceAlpha',0.7, 'MarkerEdgeAlpha',0.2);
        b.EdgeColor = group_color{gl};       
        b.LineWidth = 5;
        errorbar(gl-0.5, group_means, group_sems, 'color', group_color{gl}, 'LineWidth',3, 'CapSize',20, 'LineStyle','none');
    end
    box off
    set(gca, 'FontSize', 15, 'FontWeight', 'bold')
    set(gca,'LineWidth', 2)
    ylim([0 1])
    xticks(0.5:3);
    xticklabels(group_name);

    % show the beta of different models
    disp(' ')
    group_name = {'VTC-dependent', 'random distribution'};
    model_list = {'Category', 'TargetReward', 'Context'};
    for mi = 1:3
        mix_model{mi} = [model_list{mix_lb(mi, 1)}, '+', model_list{mix_lb(mi, 2)}];
    end
    mix_model{4} = [model_list{1}, '+', model_list{2}, '+', model_list{3}];
    model_list = [model_list, mix_model];
    for gl = 1:2
        disp('signrank: model beta>0 ')
        group_beta = [model_beta(group_label==gl, :)];
        [p,h,stats] = signrank(group_beta(:, 2), 0, 'tail', 'right');
        for mi = 1:length(model_list)
            disp([group_name{gl}, ' ', model_list{mi}, ' p-value:', num2str(p)])
        end
        disp(' ')
        disp(['figure : ', group_name{gl} ' model beta'])
        disp(' ')
        figure;
        for mi = 1:length(model_list)
            gm_beta = group_beta(:, mi);
            gm_means = mean(gm_beta);
            gm_sems  = std(gm_beta)/sqrt(size(gm_beta, 1));
            b = bar(mi-0.5, gm_means, 0.6, 'FaceColor','none', 'EdgeColor','none');
            hold on
            b.EdgeColor = group_color{gl};
            b.LineWidth = 5;
            errorbar(mi-0.5, gm_means, gm_sems, 'color', group_color{gl}, 'LineWidth',2, 'CapSize',12, 'LineStyle','none');
        end
        box off
        set(gca, 'FontSize', 15, 'FontWeight', 'bold')
        set(gca,'LineWidth', 2)
        if gl == 1
            ylim([0 0.6])
        else
            ylim([0 0.3])
        end
        xticks(0.5:7);
        xticklabels(model_list);
    end


    % Average object relations maxtrix
    disp(' ')
    disp('Group average object relations maxtrix')
    group_name = {'VTC_depended_group', 'Random_distribution_group'};
    for gl = 1:2
        disp(['figure : ', group_name{gl} ' average object relations maxtrix'])
        group_average_structure = squeeze(mean(2-estimation_structure(group_label==gl, :, :)));
        aps = (group_average_structure)./2;
        [outputArg1] = s001_showRDM(aps, it4);
    end

    % Average object value maxtrix
    disp(' ')
    disp('Group average object excepted value maxtrix')
    group_name = {'VTC_depended_group', 'Random_distribution_group'};
    for gl = 1:2
        disp(['figure : ', group_name{gl} ' average object excepted value maxtrix'])
        group_average_structure = squeeze(mean(2-rearranged_structure(group_label==gl, :, :)));
        aps = (group_average_structure)./2;
        [outputArg1] = s001_showRDM(aps, it4);
    end



    % explain the day 03 task behavior through the strcture model
    [sub_task_choice_mat] = s004_sub_task_choice_mat(Re_analysis_data_path,exp_label);
    reward_list = unique(sub_task_choice_mat{1, 1}(:, 2));
    sub_rel_value_hit = nan(size(sub_task_choice_mat, 1), 6);
    sub_irrel_value_hit = nan(size(sub_task_choice_mat, 1), 6);
    for sub_i = 1:size(sub_task_choice_mat, 1)
        sub_choice_mat = sub_task_choice_mat{sub_i, 1};
        for rel_value = 1:6
            ctxA_hit = sub_choice_mat(sub_choice_mat(:, 2)==reward_list(rel_value),  6);
            ctxB_hit = sub_choice_mat(sub_choice_mat(:, 3)==reward_list(rel_value),  7);
            sub_rel_value_hit(sub_i, rel_value) = mean([ctxA_hit; ctxB_hit]);
        end

        for irrel_value = 1:6
            ctxA_hit = sub_choice_mat(sub_choice_mat(:, 3)==reward_list(irrel_value),  6);
            ctxB_hit = sub_choice_mat(sub_choice_mat(:, 2)==reward_list(irrel_value),  7);
            sub_irrel_value_hit(sub_i, irrel_value) = mean([ctxA_hit; ctxB_hit]);
        end
    end

    STRUCT_flag = true;
    [sub_structure_model_objectEV] = s005_expected_value_of_object(Re_analysis_data_path, exp_label, STRUCT_flag, all_behavior_model);

    for sub_i = 1:size(sub_task_choice_mat, 1)
        model_theta = all_behavior_model{sub_i, 2};
        hit_p = 1./(1+exp(- model_theta(1, 4).*(sub_structure_model_objectEV{sub_i, 1}(end, :))));

        trial_reward = sub_traget_reward{sub_i, 1};
        model_choice_mat = [trial_reward, [hit_p(1:12)', hit_p(13:24)']];
        for rel_value = 1:6
            ctxA_hit = model_choice_mat(model_choice_mat(:, 1)==reward_list(rel_value),  3);
            ctxB_hit = model_choice_mat(model_choice_mat(:, 2)==reward_list(rel_value),  4);
            model_rel_value_hit(sub_i, rel_value) = mean([ctxA_hit; ctxB_hit]);
        end
    
        for irrel_value = 1:6
            ctxA_hit = model_choice_mat(model_choice_mat(:, 2)==reward_list(irrel_value),  3);
            ctxB_hit = model_choice_mat(model_choice_mat(:, 1)==reward_list(irrel_value),  4);
            model_irrel_value_hit(sub_i, irrel_value) = mean([ctxA_hit; ctxB_hit]);
        end
    end

    group_color = {[255, 128, 0]./256, [7, 126, 151]./256};
    model_color = 0.7;
    LineWidth = 2.5;
    MarkerSize = 8;
    sub_model_corr = nan(2, 2);
    disp(' ')
    disp('Figure : explain the day 03 task behavior through the strcture model')
    for gl = 1:2
        disp(['Figure :', ':', group_name{gl}])
        f = figure;
        f.Position(3:4) = [700 500];
        mean_rel_value = mean(model_rel_value_hit(group_label==gl, :));
        rel_sigmas = s011_fitsigmoid([reward_list, mean_rel_value']);
        x = -9.7:0.1:9.7;
        y = rel_sigmas(1, 3) + (1-rel_sigmas(1, 3).*2) ./ (rel_sigmas(1, 2) + exp(-rel_sigmas(1, 1) * (x)));
        plot(x, y, 'LineWidth',LineWidth, 'Color', ones(1, 3).*model_color)
        hold on
        sem=std(model_rel_value_hit(group_label==gl, :))./sqrt(sum(group_label==gl));
        errorbar(reward_list,mean_rel_value',sem, "o", 'MarkerSize',MarkerSize, 'CapSize',10, 'LineWidth',LineWidth, 'Color', ones(1, 3).*model_color)
        hold on

        mean_irrel_value = mean(model_irrel_value_hit(group_label==gl, :));
        irrel_sigmas = s011_fitsigmoid([reward_list, mean_irrel_value']);
        y = irrel_sigmas(1, 3) + (1-irrel_sigmas(1, 3).*2) ./ (irrel_sigmas(1, 2) + exp(-irrel_sigmas(1, 1) * (x)));
        plot(x, y,'--', 'LineWidth',LineWidth, 'Color', ones(1, 3).*model_color)
        hold on
        sem=std(model_irrel_value_hit(group_label==gl, :))./sqrt(sum(group_label==gl));
        errorbar(reward_list,mean_irrel_value',sem, "o", 'MarkerSize',MarkerSize, 'CapSize',10, 'LineWidth',LineWidth, 'Color', ones(1, 3).*model_color)
        hold on

        mean_rel_value = mean(sub_rel_value_hit(group_label==gl, :));
        rel_sigmas = s011_fitsigmoid([reward_list, mean_rel_value']);
        x = -9.7:0.1:9.7;
        y = rel_sigmas(1, 3) + (1-rel_sigmas(1, 3).*2) ./ (rel_sigmas(1, 2) + exp(-rel_sigmas(1, 1) * (x)));
        plot(x, y, 'LineWidth',LineWidth, 'Color', group_color{gl})
        hold on
        sem=std(sub_rel_value_hit(group_label==gl, :))./sqrt(sum(group_label==gl));
        errorbar(reward_list,mean_rel_value',sem, "o", 'MarkerSize',MarkerSize, 'CapSize',10, 'LineWidth',LineWidth, 'Color', group_color{gl})
        hold on

        mean_irrel_value = mean(sub_irrel_value_hit(group_label==gl, :));
        irrel_sigmas = s011_fitsigmoid([reward_list, mean_irrel_value']);
        y = irrel_sigmas(1, 3) + (1-irrel_sigmas(1, 3).*2) ./ (irrel_sigmas(1, 2) + exp(-irrel_sigmas(1, 1) * (x)));
        plot(x, y,'--', 'LineWidth',LineWidth, 'Color', group_color{gl})
        hold on
        sem=std(sub_irrel_value_hit(group_label==gl, :))./sqrt(sum(group_label==gl));
        errorbar(reward_list,mean_irrel_value',sem, "o", 'MarkerSize',MarkerSize, 'CapSize',10, 'LineWidth',LineWidth, 'Color', group_color{gl})
        hold on

        ylim([0 1])
        xlim([-11 11])
        reward_lim = [-9.7, -7.1, -2.6, 0, 2.6, 7.1, 9.7];
        xticks(reward_lim)
        set(gca, 'FontSize', 20, 'FontWeight', 'bold')
        set(gca,'LineWidth', LineWidth)
        box off

        sub_model_corr(gl, 1) = corr(reshape(sub_rel_value_hit(group_label==gl, :), [], 1), reshape(model_rel_value_hit(group_label==gl, :), [], 1));
        sub_model_corr(gl, 2) = corr(reshape(sub_irrel_value_hit(group_label==gl, :), [], 1), reshape(model_irrel_value_hit(group_label==gl, :), [], 1));
    end

    sub_s = [];
    for sub_i = 1:size(sub_task_choice_mat, 1)
        sub_rel_value = sub_rel_value_hit(sub_i, :);
        rel_sigmas = s006_fitsigmoid([reward_list, sub_rel_value']);
        sub_s(sub_i, 1) = rel_sigmas(1, 1);

        sub_irrel_value = sub_irrel_value_hit(sub_i, :);
        irrel_sigmas = s006_fitsigmoid([reward_list, sub_irrel_value']);
        sub_s(sub_i, 2) = irrel_sigmas(1, 1);
    end
    disp(' ')
    disp('Fitting behavioral data with sigmoid function')
    disp('rel_sigmas > irrel_sigmas, signrank')
    for gl = 1:2
        [p,h,stats] = signrank(sub_s(group_label==gl, 1), sub_s(group_label==gl, 2), 'tail', 'right');
        disp([group_name{gl}, ' p-value:', num2str(p)])
    end
    


    % Adding different prior structures to the structural model
    disp(' ')
    disp('Adding different prior structures to the structural model')
    cc_name = {'grid', 'orth', 'paral', 'VTC f1', 'VTC f2'};
    cc_angle = (15:30:360)';
    cc_coord = [cosd(cc_angle), sind(cc_angle)];
    cc_grid = [cc_coord; cc_coord];
    cc_structure(1, :) = 1-normalize(pdist(cc_grid), 'range').*2;
    cc_orth = [[cc_coord(:, 1), zeros(12, 1)]; [zeros(12, 1), cc_coord(:, 2)]];
    cc_structure(2, :) = 1-normalize(pdist(cc_orth), 'range').*2;
    ra = -90;
    rotation = [[cosd(ra); sind(ra)], [-sind(ra); cosd(ra)]];
    cc_para = [[cc_coord(:, 1), zeros(12, 1)]; [zeros(12, 1), cc_coord(:, 2)]*rotation];
    cc_structure(3, :) = 1-normalize(pdist(cc_para), 'range').*2;
    cc_f1 = [[cc_coord(:, 1), zeros(12, 1)]; [cc_coord(:, 1), zeros(12, 1)]];
    cc_structure(4, :) = 1-normalize(pdist(cc_f1), 'range').*2;
    cc_f2 = [[zeros(12, 1), cc_coord(:, 2)]; [zeros(12, 1), cc_coord(:, 2)]];
    cc_structure(5, :) = 1-normalize(pdist(cc_f2), 'range').*2;
    cc_model_NeglogLik = [];
    for cc_i = 1:5
        rep_num = 256;
        lr_LB = 0.01;
        STRUCT_flag = false;
        decayed_lr_flag = true;
        for sub_i = 1:36
            tic
            if cc_i == 3 && rotation_direction(sub_i, 1) == 1
                ra = 90;
                rotation = [[cosd(ra); sind(ra)], [-sind(ra); cosd(ra)]];
                cc_para = [[cc_coord(:, 1), zeros(12, 1)]; [zeros(12, 1), cc_coord(:, 2)]*rotation];
                prior_structure = 1-normalize(pdist(cc_para), 'range').*2;
            else
                prior_structure = cc_structure(cc_i, :);
            end

            sub_choices = cell2mat(sub_learning_trials{sub_i, 2}');
            stimuli_id = (sub_choices(:, 2)-1)*12+sub_choices(:, 3);
            outcomes = sub_choices(:, 4);
            choices = sub_choices(:, 5);

            if decayed_lr_flag
                % [learning_rate, decay_rate, decay_steps, beta]
                LB = zeros(1, 4); % Lower bound
                UB = [1, 1, 1080, 10.*ones(1, 1)]; % Upper bound
                X0 = LB + (UB-LB).*rand(rep_num, size(LB, 2)); % random initialisation
            else
                % [learning_rate, beta]
                LB = zeros(1, 2); % Lower bound
                UB = [lr_LB, 10.*ones(1, 1)]; % Upper bound
                X0 = LB + (UB-LB).*rand(rep_num, size(LB, 2)); % random initialisation
            end
            energy = @s003_ctx_RL_LogLikelihood;
            OPTIM_options = optimset('Display', 'off') ;

            NeglogLik_ML_all = nan(rep_num,1);
            theta_ML_all= nan(rep_num, size(LB, 2)); % [alpha, beta, stimuli_structure]
            parfor ri = 1:rep_num
                [theta_ML_all(ri,:), NeglogLik_ML_all(ri, 1)] = fmincon(@(theta) energy(theta,stimuli_id,outcomes,choices,STRUCT_flag, decayed_lr_flag, prior_structure),X0(ri,:)',[],[],[],[],LB', UB',[],OPTIM_options);
            end
            theta_all = [theta_ML_all, NeglogLik_ML_all];
            theta_ML = sortrows(theta_all, size(theta_all, 2));
            X0 = sortrows([X0, NeglogLik_ML_all], size(theta_all, 2));
            cc_model_NeglogLik(sub_i, cc_i) = theta_ML(1, end);

            t = toc;
            disp([cc_name{cc_i} '_', sub_dir(sub_i).name, '_t', num2str(t)])
        end
    end

     for gl = 1:2
        disp(' ')
        disp(['Figure :', group_name{gl}])
        group_NeglogLik = cc_model_NeglogLik(group_label==gl, 1)-cc_model_NeglogLik(group_label==gl, 2:end);
        figure;
        for mi = 1:size(group_NeglogLik, 2)
            gm_NeglogLik = group_NeglogLik(:, mi);
            gm_means = mean(gm_NeglogLik);
            gm_sems  = std(gm_NeglogLik)/sqrt(size(gm_NeglogLik, 1));
            jitter = 0.2;
            b = bar(mi-0.5, gm_means, 0.6, 'FaceColor','none', 'EdgeColor','none');
            hold on
            xi = mi-0.5 + (rand(size(gm_NeglogLik)) - 0.5) * 2 * jitter;
            scatter(xi, gm_NeglogLik, 50, ones(1, 3).*0.2, 'filled', ...
                'MarkerFaceAlpha',0.7, 'MarkerEdgeAlpha',0.2);
            b.EdgeColor = group_color{gl};
            b.LineWidth = 5;
            errorbar(mi-0.5, gm_means, gm_sems, 'color', group_color{gl}, 'LineWidth',2, 'CapSize',12, 'LineStyle','none');
        end
        box off
        set(gca, 'FontSize', 15, 'FontWeight', 'bold')
        set(gca,'LineWidth', 2)
        if gl == 1
            ylim([-400 400])
        else
            ylim([-150 100])
        end
        xticks(0.5:7);
        xticklabels(cc_name(2:end));

        disp([group_name{gl}, ', signrank:'])
        for mi = 1:4
            [p,h,stats] = signrank(group_NeglogLik(:, mi), 0, 'tail','right');
            disp([cc_name{mi+1}, '>grid, p-value:', num2str(p)])
        end
    end

    outputArg1 = 'Analysis completed without error';
end

