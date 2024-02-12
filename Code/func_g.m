function [out] = func_g(t,x)
%FUNC_G Summary of this function goes here
%   Detailed explanation goes here
out=(1/2+1/2*cos(pi*x))*(abs(x)<=1)*(0<=t)*(t<=1);
end

