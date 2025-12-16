function [R] = s010_circular_gaussian_function(x, theta)
%s010_circular_gaussian_function 此处显示有关此函数的摘要
%   此处显示详细说明
% von Mises distribution
        R0 = x(1); % Baseline response
        A = x(2); % Amplitude
        mu = x(3); % Preferred direction
        sigma = x(4); % Tuning width
        R = R0 + A * exp((cos(theta - mu) - 1) / (2 * sigma^2));
end



% 
% theta = -pi:0.1:pi;
% R0 = 5; % Baseline response
% A = 3; % Amplitude
% mu = -pi/3; % Preferred direction
% sigma = 0.5; % Tuning width
% R = R0 + A * exp((cos(theta - mu) - 1) / (2 * sigma^2));
% figure
% polarplot(theta,R)
% 








