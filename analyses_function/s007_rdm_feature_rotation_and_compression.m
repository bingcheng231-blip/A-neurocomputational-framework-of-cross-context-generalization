function [roi_features_dist, roi_rotation_rr] = s007_rdm_feature_rotation_and_compression(roi_cc_rdm, roi_num, sub_feature_list)
%S007_RDM_FEATURE_ROTATION_AND_COMPRESSION 此处显示有关此函数的摘要
%   此处显示详细说明

    rep_num = 32;
    LB = [0, zeros(1, 4)]; % Lower bound
    UB = [10, ones(1, 4)]; % Upper bound
    params_num = length(LB);
    OPTIM_options = optimset('Display', 'off') ;
    for sub_i = 1:size(roi_cc_rdm, 1)
        tic
        % Representational geometry
        target_rdm = zscore(roi_cc_rdm{sub_i, roi_num, 4}');
        sub_feature = sub_feature_list{sub_i, 1};
        parfor ra_i = 1:361
            rotation = [[cosd(ra_i); sind(ra_i)], [-sind(ra_i); cosd(ra_i)]];
            X0 = LB + (UB-LB).*rand(rep_num, params_num); % random initialisation
            theta_ML_all = nan(rep_num, params_num); NeglogLik_ML_all = nan(rep_num, 1);
            for ri = 1:rep_num
                lastwarn('');
                [theta_ML_all(ri,:), NeglogLik_ML_all(ri)] = fmincon(@(theta) sumsqr(target_rdm'-s008_ROIs_rdm_ctxDist_calc(sub_feature, rotation, theta)),X0(ri,:),[],[],[],[],LB', UB',[],OPTIM_options);
                [warnMsg, warnId] = lastwarn;
                rep_c = 0;
                while ~isempty(warnMsg)
                    disp('Recalculating...')
                    lastwarn('');
                    X0 = LB + (UB-LB).*rand(rep_num, params_num); % random initialisation
                    [theta_ML_all(ri,:), NeglogLik_ML_all(ri)] = fmincon(@(theta) sumsqr(target_rdm'-s008_ROIs_rdm_ctxDist_calc(sub_feature, rotation, theta)),X0(ri,:),[],[],[],[],LB', UB',[],OPTIM_options);
                    [warnMsg, warnId] = lastwarn;
                    rep_c = rep_c+1;
                    if rep_c>100
                        error('Error initialization')
                    end
                end
            end
            theta_all = [theta_ML_all, NeglogLik_ML_all];
            theta_ML = sortrows(theta_all, size(theta_all, 2));

            roi_features_dist(sub_i, ra_i, :) = theta_ML(1, 1:end-1);
            roi_rotation_rr(sub_i, ra_i) = corr(target_rdm, s008_ROIs_rdm_ctxDist_calc(sub_feature, rotation, theta_ML(1, 1:end-1))');
        end
    end
    t = toc;
    disp(['Sub_', num2str(sub_i, '%02d'), '_t', num2str(t)])

end

