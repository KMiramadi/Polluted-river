function [c_j] = func_c_j(eta,h,x_val)
%FUNC_C_J Summary of this function goes here
%   Detailed explanation goes here

c_j=eta/h^2+1/(func_w(x_val)*2*h);
end

