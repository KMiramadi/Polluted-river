function [a_j] = func_a_j(eta,h,x_val)
%FUNC_A_J Summary of this function goes here
%   Detailed explanation goes here
a_j=eta/h^2-1/(func_w(x_val)*2*h);

end

