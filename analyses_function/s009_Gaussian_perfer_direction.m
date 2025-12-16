function [x_list] = s009_Gaussian_perfer_direction(responses, v_theta0, theta, rep_num)
%s009_Gaussian_perfer_direction 此处显示有关此函数的摘要
%   此处显示详细说明
% Initial parameters and bounds

   % x0 = [Baseline response; Amplitude; Preferred direction; Tuning width]
   if isempty(v_theta0)
       theta0 = -pi+2*rand(rep_num)*pi;
   else
       theta0 = v_theta0+(rand(rep_num)-0.5)*2*pi;
   end

    
    lb = [0, 0, -pi, 0];
    ub = [inf, inf, pi, inf];
    
    % Fit model
    options = optimset('Display','off');
    x_list = nan(rep_num, 5);
    parfor rep_i = 1:rep_num
        x0 = [mean(responses), max(responses) - min(responses), theta0(rep_i), 1];
        [x, resnorm] = lsqcurvefit(@s010_circular_gaussian_function, x0, theta, responses, lb, ub, options);
        x_list(rep_i, :) = [x, resnorm];
    end
    x_list = sortrows(x_list, 5);


end











