function [output] = func_w(x)
%FUNC_W Summary of this function goes here
%   Detailed explanation goes here

output=1+0.5.*exp(-100.*(x-0.3).^2)-4/5.*exp(-250*(x-0.6).^2)+(1-x)./5;

end

