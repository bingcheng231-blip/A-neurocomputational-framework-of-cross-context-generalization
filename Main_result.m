clear all
it4=wd_compit4_1;
close all



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% setting before running  %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make sure you have already added the function folder ./analyses_function
% and the utilities folder ./analyses_function/utilities to your MATLAB path.

% chang this path to yuor ./data path;
Re_analysis_data_path = 'D:\Context_learning\demo\data\';

% chang this path to yuor ./figure_results;
figure_generation_data_path = 'D:\Context_learning\demo\figure_results\'; 

% swith to ture to re run all main analysis
re_analysis = false; 


%%%%%%%%%%%%%%%%%%%%%%%%%    Scripts info    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scripts for A neurocomputational framework of cross-context generalization: dynamic 
% representational geometry in high-level visual cortex and dual coding in vmPFC
%
% First Authors: Bincheng Wen
% Correspondence Authors: Le Chang (lechang@ion.ac.cn); Huiguang He (huiguang.he@ia.ac.cn)
%
% In general, if anything across the codes is missing or unlclear, you are
% most welcome to contact Huiguang He, always happy to chat and explain whatever 



tic
%% 1.1 stimulus distribution (Supplemental figure 1)
% close all

if re_analysis
    [outputArg1] = f001_stimulus_distribution(Re_analysis_data_path);
    disp(outputArg1)
else
    load(fullfile(figure_generation_data_path, 'r001_stimulus_rdm.mat'))
    stimulus_PCscore = load([Re_analysis_data_path, 'd001_Normalized_stimuluse_PCs.mat']);
    Experiment_stimuli = load([Re_analysis_data_path, 'd002_Stimulus_image_numbers_for_each_category.mat']).Experiment_stimuli;
    
    rdm_name = {'ideal circular', 'all experiment stimuli', 'day 02 learning stimuli',...
        'day 01 localizer stimuli', 'day 03 task stimuli'};

    % Correlation between the ideal circular and the experiment stimuli
    disp('Figure : Correlation between the ideal circular and the experiment stimuli')
    mdl = fitlm(zscore(all_line_rdm(1, :)'), zscore(all_line_rdm(2, :)'));
    [rho,pval] = corr(all_line_rdm(1, :)', all_line_rdm(2, :)', 'Type', 'Spearman');
    disp(['Correlation:', num2str(rho), ' p:', num2str(pval)])
    f = figure;
    f.Position(3:4) = [500 500];
    point_color = [0.5, 0.5, 0.5];
    scatter(zscore(all_line_rdm(1, :)'), zscore(all_line_rdm(2, :)'), 35, point_color, 'filled');
    hold on
    x_range = -2:0.1:1.5;
    predict_line = predict(mdl, x_range')';
    plot(x_range, predict_line,'Color', [1 0 0], 'linewidth', 3)
    set(gca, 'FontSize', 25, 'FontWeight', 'bold')
    set(gca,'LineWidth',3)

    % RDM for the 12 centers of each stimulus set
    for rdm_i = 1:5
        showRDM = squareform(all_line_rdm(rdm_i, :));
        aps = showRDM./max(showRDM, [], 'all');
        disp(['Figure : Show the RDM for the 12 centers of ', rdm_name{rdm_i}])
        [outputArg1] = s001_showRDM(aps, it4);
    end

    % cross stimulus_sets correlation
    disp('Figure : cross stimulus_sets correlation')
    set_corr_rdm = corr(all_line_rdm', 'Type', 'Spearman')-0.8;
    aps = set_corr_rdm./max(set_corr_rdm, [], 'all');
    [outputArg1] = s001_showRDM(aps, it4);


    % show example stimulus for each category
    im_show = [];
    for ci = 1:12
        center_idex = strcmp(stimulus_PCscore.img_list(:, 1), ['center_', num2str(ci)]);
        center_stimuli_list = stimulus_PCscore.img_list(center_idex, :);
        im_show{1, ci} = imresize(center_stimuli_list{1, 3}, [224, 224]);
    end
    disp('Figure : example stimulus for each category')
    figure
    imshow(cell2mat(im_show))
end




%% 1.2 participant learning behavior analysis (figure 1)
% close all


if re_analysis
    % Average acc of each run during the learning
    exp_label = 'all'; % exp_label = ['all', 'MEG', 'fMRI'], for participants
    % who participated in all experiments, only MEG experiments, or only fMRI 
    % experiments, respectively.
    [outputArg1] = f002_sub_learning_acc(Re_analysis_data_path, exp_label);

    % Average acc of the task
    [outputArg1] = f004_sub_task_acc(Re_analysis_data_path, exp_label);

    % Modeling learning behavior 
    [outputArg1] = f003_behavior_learning_model(Re_analysis_data_path, exp_label);

    % Modeling learning behavior
    [outputArg1] = f005_task_choice_RSA(Re_analysis_data_path, exp_label);
    disp(outputArg1)


else
    % Average acc of each run during the learning
    disp(' ')
    load(fullfile(figure_generation_data_path, 'r002_learning_session_acc.mat'))
    exp_label = 'all'; % exp_label = ['all', 'MEG', 'fMRI'], for participants
    % who participated in all experiments, only MEG experiments, or only fMRI 
    % experiments, respectively.
    [group_label, ~, sub_traget_reward, ~] = s002_grouplabel(Re_analysis_data_path, exp_label);

    % The original figure was generated using GraphPad.
    % VTC_depended_group
    g01_lr_acc = learning_session(:, group_label==1);
    % Random_distribution_group
    g02_lr_acc = learning_session(:, group_label==2);
    disp('Figure : Average acc of each run during the learning experiment')
    figure
    group_color = {[255, 128, 0]./256, [7, 126, 151]./256};
    y = g01_lr_acc'; x = 1:15;
    SEM=std(y)./sqrt(size(y, 1));
    s = shadedErrorBar(x,mean(y, 1),SEM, 'lineProps', '-');
    s.patch.FaceColor = [0.5, 0.5, 0.5];
    set(s.edge,'LineWidth',0.5, 'Color', [0.7, 0.7, 0.7])
    hold on
    plot(x, mean(y, 1), 'LineWidth', 2, 'color', group_color{1})
    hold on
    y = g02_lr_acc'; x = 1:15;
    SEM=std(y)./sqrt(size(y, 1));
    s = shadedErrorBar(x,mean(y, 1),SEM, 'lineProps', '-');
    s.patch.FaceColor = [0.5, 0.5, 0.5];
    set(s.edge,'LineWidth',0.5, 'Color', [0.7, 0.7, 0.7])
    hold on
    plot(x, mean(y, 1), 'LineWidth', 2, 'color', group_color{2})
    box off
    set(gca, 'FontSize', 15, 'FontWeight', 'bold')
    set(gca,'LineWidth', 2)

    % Average acc of the task
    disp(' ')
    load(fullfile(figure_generation_data_path, 'r003_task_session_acc.mat'))
    % VTC_dependent_group
    group_acc{1} = mean(task_session(:, group_label==1));
    % Random_distribution_group
    group_acc{2} = mean(task_session(:, group_label==2));
    [p,h,stats] = ranksum(group_acc{1}, group_acc{2}, 'tail', 'right');
    disp(['task acc, ranksum, VTC dependent > Random distribution: pvalue=', num2str(p)])
    disp('Figure : Average acc of the task experiment')
    group_name = {'VTC-dependent', 'random distribution'};
    group_color = {[255, 128, 0]./256, [7, 126, 151]./256};
    figure;
    for gl = 1:2
        group_data = group_acc{gl};
        group_means = mean(group_data);
        group_sems  = std(group_data)/sqrt(size(group_data, 1));
        jitter = 0.2;
        b = bar(gl-0.5, group_means, 0.6, 'FaceColor','none', 'EdgeColor','none');
        hold on
        xi = gl-0.5 + (rand(size(group_data)) - 0.5) * 2 * jitter;
        scatter(xi, group_data, 50, ones(1, 3).*0.2, 'filled', ...
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


    % behavior model comparison
    disp(' ')
    load(fullfile(figure_generation_data_path, 'r004_behavior_model.mat'))
    model_NeglogLik = [];
    for sub_i = 1:size(all_behavior_model, 1)
        model_theta = all_behavior_model{sub_i, 1};
        model_NeglogLik(sub_i, 1) = model_theta(1, end);

        model_theta = all_behavior_model{sub_i, 2};
        model_NeglogLik(sub_i, 2) = model_theta(1, end);
    end
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
    % show the object relations repeatability
    disp(' ')
    disp('figure : the repeatability of the object relations derived from the structure model')
    load(fullfile(figure_generation_data_path, 'r006_object_relations_from_the_structure_model.mat'))
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
    mix_lb = nchoosek(1:3, 2);
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


    load(fullfile(figure_generation_data_path, 'r007_sub_task_sigmoid.mat'))
    disp(' ')
    disp('Fitting behavioral data with sigmoid function')
    disp('rel_sigmas > irrel_sigmas, signrank')
    for gl = 1:2
        [p,h,stats] = signrank(sub_s(group_label==gl, 1), sub_s(group_label==gl, 2), 'tail', 'right');
        disp([group_name{gl}, ' p-value:', num2str(p)])
    end


    % Adding different prior structures to the structural model
    disp(' ')
    load(fullfile(figure_generation_data_path, 'r005_behavior_model_with_different_prior_structures.mat'))
    cc_name = {'grid', 'orth', 'paral', 'VTC f1', 'VTC f2'};
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


    % task choice RSA
    disp('figure : task choice RSA ')
    load(fullfile(figure_generation_data_path, 'r008_task_choice_RSA.mat'))
    disp(' ')
    figure;
    group_color = {[255, 128, 0]./256, [7, 126, 151]./256};
    group_beta = sub_beta(group_label==1, :);
    figure;
    gl = 1;
    for mi = 1:size(group_beta, 2)
        gm_beta = group_beta(:, mi);
        gm_means = mean(gm_beta);
        gm_sems  = std(gm_beta)/sqrt(size(gm_beta, 1));
        jitter = 0.2;
        b = bar(mi-0.5, gm_means, 0.6, 'FaceColor','none', 'EdgeColor','none');
        hold on
        xi = mi-0.5 + (rand(size(gm_beta)) - 0.5) * 2 * jitter;
        scatter(xi, gm_beta, 50, ones(1, 3).*0.2, 'filled', ...
            'MarkerFaceAlpha',0.7, 'MarkerEdgeAlpha',0.2);
        b.EdgeColor = group_color{gl};
        b.LineWidth = 5;
        errorbar(mi-0.5, gm_means, gm_sems, 'color', group_color{gl}, 'LineWidth',2, 'CapSize',12, 'LineStyle','none');
    end
    gl = 2;
    group_beta = sub_beta(group_label==gl, :);
    for mi = 1:size(group_beta, 2)
        gm_beta = group_beta(:, mi);
        gm_means = mean(gm_beta);
        gm_sems  = std(gm_beta)/sqrt(size(gm_beta, 1));
        jitter = 0.2;
        b = bar(mi-0.5+2, gm_means, 0.6, 'FaceColor','none', 'EdgeColor','none');
        hold on
        xi = mi-0.5+2 + (rand(size(gm_beta)) - 0.5) * 2 * jitter;
        scatter(xi, gm_beta, 50, ones(1, 3).*0.2, 'filled', ...
            'MarkerFaceAlpha',0.7, 'MarkerEdgeAlpha',0.2);
        b.EdgeColor = group_color{gl};
        b.LineWidth = 5;
        errorbar(mi-0.5+2, gm_means, gm_sems, 'color', group_color{gl}, 'LineWidth',2, 'CapSize',12, 'LineStyle','none');
    end
    box off
    set(gca, 'FontSize', 15, 'FontWeight', 'bold')
    set(gca,'LineWidth', 2)
    ylim([-0.5 1])
    xticks(0.5:4);
    xticklabels({'Boundary', 'Linear', 'Boundary', 'Linear'});
end







%% 2.1. Dual representation of value and category in the vmPFC (figure 2)
% close all

if re_analysis
    % Average acc of each run during the learning
    exp_label = 'fMRI'; % exp_label = ['all', 'MEG', 'fMRI'], for participants
    % who participated in all experiments, only MEG experiments, or only fMRI 
    % experiments, respectively.

    % vmPFC model based RSA
    [outputArg1] = f006_vmPFC_model_RDM_fitting(Re_analysis_data_path,exp_label);

    % vmPFC representational geometry
    [outputArg1] = f007_vmPFC_feature_rotation_and_compression(Re_analysis_data_path,exp_label);


else
    % Average acc of each run during the learning
    exp_label = 'fMRI'; % exp_label = ['all', 'MEG', 'fMRI'], for participants
    % who participated in all experiments, only MEG experiments, or only fMRI 
    % experiments, respectively.
    [group_label, ~, sub_traget_reward, ~] = s002_grouplabel(Re_analysis_data_path, exp_label);

    % vmPFC model based RSA
    [outputArg1] = f006_vmPFC_model_RDM_fitting(Re_analysis_data_path,exp_label);

    % vmPFC representational geometry
    disp(' ')
    disp('figure : Geometric relationships of category’s EVs across different contexts in the vmPFC on Day 3')
    load(fullfile(figure_generation_data_path, 'r009_fMRI_ROIs_EV_representational_geometry.mat'))
    roi_id = 2; % vmPFC
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
        rotation_rr = squeeze(roi_vv_rotation_rr(group_label==gl, roi_id, :));
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
    load(fullfile(figure_generation_data_path,  'r010_vmPFC_mds_rdm.mat'))
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


    disp(' ')
    load(fullfile(figure_generation_data_path,  'r011_fMRI_ROIs_theoretical_EV_representation.mat'))
    disp('Group level cross-context EVs rotation in the vmPFC:')
    roi_num = 2;
    for gl = 1:2
        disp(group_name{gl})
        group_ra = sub_roi_ra(group_label==gl, roi_num);
        mean_rotation = mean(group_ra);
        SEM_rotation=std(group_ra)./sqrt(size(group_ra, 1));
        disp(['mean_rotation: ', num2str(mean_rotation)]);
        disp(['sem: ', num2str(SEM_rotation)])
        disp(' ')
    end


    load(fullfile(figure_generation_data_path,  'r010_fMRI_ROIs_VTCfeatures_representational_geometry.mat'))
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
        rotation_rr = squeeze(roi_cc_rotation_rr(group_label==gl, roi_id, :));
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


    disp(' ')
    load(fullfile(figure_generation_data_path,  'r012_fMRI_ROIs_theoretical_VTCfeatures_representation.mat'))
    disp('Group level cross-context VTC features rotation in the vmPFC:')
    roi_num = 2;
    for gl = 1:2
        disp(group_name{gl})
        group_ra = sub_roi_ra(group_label==gl, roi_num);
        mean_rotation = mean(group_ra);
        SEM_rotation=std(group_ra)./sqrt(size(group_ra, 1));
        disp(['mean_rotation: ', num2str(mean_rotation)]);
        disp(['sem: ', num2str(SEM_rotation)])
        disp(' ')
    end

end




%% 3.1. Phase and amplitude modulation in VTC (figure 3)
% close all


% Average acc of each run during the learning
exp_label = 'fMRI'; % exp_label = ['all', 'MEG', 'fMRI'], for participants
% who participated in all experiments, only MEG experiments, or only fMRI
% experiments, respectively.

% the phase and amplitude modulation of VTC voxel
[outputArg1] = f008_VTC_phase_and_amplitude_modulation(Re_analysis_data_path,exp_label);

% VTC representational geometry
[outputArg1] = f009_VTC_feature_rotation_and_compression(Re_analysis_data_path,exp_label,figure_generation_data_path, re_analysis);

  




%% 4.1. Value estimation based on VTC representations (figure 4)
% close all

% Average acc of each run during the learning
exp_label = 'all'; % exp_label = ['all', 'MEG', 'fMRI'], for participants
% who participated in all experiments, only MEG experiments, or only fMRI
% experiments, respectively.

% generate the training set for the gaussian mlp training
% the mlp training was accomplised by pytorch, see 'MLP_with_modulated_GaussianActivation_verA100.ipynb'
[outputArg1] = f010_Generate_MLP_training_data(Re_analysis_data_path,exp_label);

% loading the trained MLP results for further analysis
[outputArg1] = f011_Gaussian_MLP_simulation(Re_analysis_data_path,exp_label, re_analysis);

% ROI representational geometry
exp_label = 'fMRI'; 
[outputArg1] = f012_ROI_feature_rotation_and_compression(Re_analysis_data_path,exp_label, figure_generation_data_path, re_analysis);






%% 5.1. Value estimation based on VTC representations (figure 5)
% close all


% Average acc of each run during the learning
exp_label = 'MEG'; % exp_label = ['all', 'MEG', 'fMRI'], for participants
% who participated in all experiments, only MEG experiments, or only fMRI
% experiments, respectively.

% the temporal evolution during the stimuli presentation of the learning task
[outputArg1] = f013_MEG_sensor_level_dynamic(Re_analysis_data_path, exp_label, figure_generation_data_path);


% Dynamic interactions between different brain regions during the stimuli
% presentation of the MEG learning task.
[outputArg1] = f014_MEG_ROI_interactions(Re_analysis_data_path, exp_label, figure_generation_data_path);




t = toc;
disp(['total time: ', num2str(t)])




















