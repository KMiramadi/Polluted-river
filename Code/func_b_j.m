function [b_j] = func_b_j(eta,h,x_val)
%FUNC_B_J Summary of this function goes here
%   Detailed explanation goes here
b_j=-2*eta/h^2-func_w_derivative(x_val)/func_w(x_val)^2;
end

