function [outputArg1] = f001_stimulus_distribution(Re_analysis_data_path)
%F001_ 此处显示有关此函数的摘要
%   此处显示详细说明
    stimulus_PCscore = load([Re_analysis_data_path, 'd001_Normalized_stimuluse_PCs.mat']);
    Experiment_stimuli = load([Re_analysis_data_path, 'd002_Stimulus_image_numbers_for_each_category.mat']).Experiment_stimuli;
    it4=wd_compit4_1;
    close;

    % idealistic circular
    rn = 1;radius = 45;
    circular_center = zeros(12, 64);angel = 15;
    rlimit(rn) = radius*sind(15);
    for node_i = 1:12
        circular_center(node_i, 2) = radius.*cosd(angel);
        circular_center(node_i, 3) = radius.*sind(angel);
        angel = angel+30;
    end
    circular_rdm = pdist2(circular_center, circular_center);
    aps = circular_rdm./max(circular_rdm, [], 'all');
    disp('Figure01: Show the RDM for an ideal circular arrangement of 12 centers')
    [outputArg1] = s001_showRDM(aps, it4);
    line_circular_rdm = pdist(circular_center);
    all_line_rdm(1, :) = line_circular_rdm;
    nn_circular_rdm = zscore(line_circular_rdm);


    % stimuli circular
    activation = stimulus_PCscore.activation(2:end, :);
    center_list = cell(12,1);
    stimuli_center_list = zeros(12, 64);
    for ci = 1:12
        center_stimuli_list = cell(1, 3);
        im_n = 1;
        for im_i = 1:size(stimulus_PCscore.img_list, 1)
            im_name = strtrim(stimulus_PCscore.img_list{im_i, 2});
            im_center = str2num(im_name(end-9:end-8));
            if ci == im_center
                center_stimuli_list{im_n, 1} = stimulus_PCscore.img_list{im_i, 1};
                center_stimuli_list{im_n, 2} = im_name;
                center_stimuli_list{im_n, 3} = activation(im_i, :);
                im_n = im_n+1;
            end
        end
        stimuli_activation = cell2mat(center_stimuli_list(:, 3));
        stimuli_center_list(ci, :) = mean(stimuli_activation);
        center_list{ci, 1} = center_stimuli_list;

        stimuli_center_dist = pdist2(stimuli_center_list(ci, :), stimuli_activation);

        scd(ci, 1) = min(stimuli_center_dist);
        scd(ci, 2) = mean(stimuli_center_dist);
        scd(ci, 3) = max(stimuli_center_dist);
    end
    stimuli_center_rdm = pdist2(stimuli_center_list, stimuli_center_list);
    aps = stimuli_center_rdm./max(stimuli_center_rdm, [], 'all');
    disp('Figure02: Show the RDM for the 12 centers of all the experiment stimuli')
    [outputArg1] = s001_showRDM(aps, it4);
    line_stimuli_center_rdm = pdist(stimuli_center_list);
    all_line_rdm(2, :) = line_stimuli_center_rdm;
    nn_stimuli_center_rdm = zscore(line_stimuli_center_rdm);
    
    disp('Figure03: Correlation between the ideal circular and the experiment stimuli')
    mdl = fitlm(nn_circular_rdm, nn_stimuli_center_rdm);
    [rho,pval] = corr(nn_circular_rdm', nn_stimuli_center_rdm', 'Type', 'Spearman');
    disp(['Correlation:', num2str(rho), ' p:', num2str(pval)])
    f = figure;
    f.Position(3:4) = [500 500];
    point_color = [0.5, 0.5, 0.5];
    scatter(nn_circular_rdm', nn_stimuli_center_rdm', 35, point_color, 'filled');
    hold on
    x_range = -2:0.1:1.5;
    predict_line = predict(mdl, x_range')';
    plot(x_range, predict_line,'Color', [1 0 0], 'linewidth', 3)
    set(gca, 'FontSize', 25, 'FontWeight', 'bold')
    set(gca,'LineWidth',3)

    %localizer set distribution
    localizer_center = zeros(12, 64);dist_list = zeros(12, 18);
    localizer_img_num = Experiment_stimuli{2, 2};
    rep_localizer_img_num = Experiment_stimuli{2, 3};
    for ci = 1:12
        center_stimuli_list = center_list{ci, 1};
        stimuli_activation = cell2mat(center_stimuli_list(:, 3));
        localizer_stimuli_activation = [stimuli_activation(localizer_img_num, :);stimuli_activation(rep_localizer_img_num, :)];
        localizer_center(ci, :) = mean(localizer_stimuli_activation);

        stimuli_center_dist = pdist2(stimuli_center_list(ci, :), localizer_stimuli_activation);
        dist_list(ci, :) = stimuli_center_dist;
        scd(ci, 1) = min(stimuli_center_dist);
        scd(ci, 2) = mean(stimuli_center_dist);
        scd(ci, 3) = max(stimuli_center_dist);
    end
    dist_list = dist_list';
    localizer_center_rdm = pdist2(localizer_center, localizer_center);
    aps = localizer_center_rdm./max(localizer_center_rdm, [], 'all');
    disp('Figure04: Show the RDM for the 12 centers of the day 01 localizer experiment stimuli')
    [outputArg1] = s001_showRDM(aps, it4);
    localizer_rdm = pdist(localizer_center);
    all_line_rdm(4, :) = localizer_rdm;

    
    %learning set distribution
    learning_center = zeros(12, 64); dist_list = zeros(12, 45);
    learning_img_num = Experiment_stimuli{3, 2};
    for ci = 1:12
        center_stimuli_list = center_list{ci, 1};
        stimuli_activation = cell2mat(center_stimuli_list(:, 3));
        learning_stimuli_activation = stimuli_activation(learning_img_num, :);
        learning_center(ci, :) = mean(learning_stimuli_activation);

        stimuli_center_dist = pdist2(stimuli_center_list(ci, :), learning_stimuli_activation);
        dist_list(ci, :) = stimuli_center_dist;
        scd(ci, 1) = min(stimuli_center_dist);
        scd(ci, 2) = mean(stimuli_center_dist);
        scd(ci, 3) = max(stimuli_center_dist);
    end
    dist_list = dist_list';
    learning_center_rdm = pdist2(learning_center, learning_center);
    aps = learning_center_rdm./max(learning_center_rdm, [], 'all');
    disp('Figure05: Show the RDM for the 12 centers of the day 02 learning experiment stimuli')
    [outputArg1] = s001_showRDM(aps, it4);
    learning_rdm = pdist(learning_center);
    all_line_rdm(3, :) = learning_rdm;


    %task set distribution
    task_center = zeros(12, 64);dist_list = zeros(12, 18);
    task_img_num = Experiment_stimuli{4, 2};
    rep_task_img_num = Experiment_stimuli{4, 3};
    for ci = 1:12
        center_stimuli_list = center_list{ci, 1};
        stimuli_activation = cell2mat(center_stimuli_list(:, 3));
        task_stimuli_activation = [stimuli_activation(task_img_num, :);stimuli_activation(rep_task_img_num, :)];
        task_center(ci, :) = mean(task_stimuli_activation);

        stimuli_center_dist = pdist2(stimuli_center_list(ci, :), task_stimuli_activation);
        dist_list(ci, :) = stimuli_center_dist;
        scd(ci, 1) = min(stimuli_center_dist);
        scd(ci, 2) = mean(stimuli_center_dist);
        scd(ci, 3) = max(stimuli_center_dist);
    end
    dist_list = dist_list';
    task_center_rdm = pdist2(task_center, task_center);
    aps = task_center_rdm./max(task_center_rdm, [], 'all');
    disp('Figure06: Show the RDM for the 12 centers of the day 03 task experiment stimuli')
    [outputArg1] = s001_showRDM(aps, it4);
    task_rdm = pdist(task_center);
    all_line_rdm(5, :) = task_rdm;


    disp('Figure07: cross stimulus_sets correlation')
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
    disp('Figure08: example stimulus for each category')
    figure
    imshow(cell2mat(im_show))

    outputArg1 = 'Analysis completed without error';
end

