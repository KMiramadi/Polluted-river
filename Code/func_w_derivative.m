function [derivative] = func_w_derivative(x)
%FUNC_W_DERIVATIVE Summary of this function goes here
%   Detailed explanation goes here
derivative=-100*(x-3/10)*exp(-100*(x-3/10)^2)+400*(x-3/5)*exp(-250*(x-3/5)^2)-1/5;
end

