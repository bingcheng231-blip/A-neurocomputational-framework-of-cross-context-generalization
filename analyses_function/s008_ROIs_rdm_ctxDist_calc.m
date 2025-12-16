function [roi_rdm, ctx_vxy] = s008_ROIs_rdm_ctxDist_calc(sub_vv_coord, rotation, theta)
%F082_ROIS_RDM_CTXDIST_CALC 此处显示有关此函数的摘要
%   此处显示详细说明
    % theta = X0(1, :);
    % theta(theta>1-1e-5) = 1-1e-5;
    ctx_dist = [zeros(size(sub_vv_coord, 1), 1); ones(size(sub_vv_coord, 1), 1).*theta(1)];
    ctxA_vxy = [sub_vv_coord(:, 1).*(1-theta(2)), sub_vv_coord(:, 2).*(1-theta(3))];
    ctxB_vxy = [sub_vv_coord(:, 1).*(1-theta(4)), sub_vv_coord(:, 2).*(1-theta(5))];
    ctxB_vxy = ctxB_vxy * rotation;
    ctx_vxy = [[ctxA_vxy; ctxB_vxy], ctx_dist];
    roi_rdm = zscore(pdist(ctx_vxy));
end

