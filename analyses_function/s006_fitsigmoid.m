function [theta_ML] = s006_fitsigmoid(data)
%s006_fitsigmoid 此处显示有关此函数的摘要
%   此处显示详细说明
    rep_num = 64;

    LB = zeros(1, 3); % Lower bound for [alpha, beta]
    UB = [20.*ones(1, 2), 0.5]; % Upper bound for [alpha, beta]

    X0 = [20.*rand(rep_num, 2), 0.5.*rand(rep_num, 1)]; % random initialisation

    betas_all = nan(rep_num, 3);
    NeglogLik_ML_all = nan(rep_num, 1);

    OPTIM_options = optimset('Display', 'off') ;
    
    for ri = 1:rep_num
        [betas_all(ri,:), NeglogLik_ML_all(ri)] = fmincon(@(betas)sumsqr(data(:, 2)-mysigmoid(betas, data(:, 1))),X0(ri,:)',[],[],[],[],LB', UB',[],OPTIM_options);
    end
    theta_all = [betas_all, NeglogLik_ML_all];
    theta_ML = sortrows(theta_all, size(theta_all, 2));
    function fun = mysigmoid(b, x)
        fun = b(3) + (1-b(3).*2) ./ (b(2) + exp(-b(1) * (x)));
    end

end

