function [psi] = psimaker(t,h_x)
%
xspace=0:h_x:1;
psi=zeros(length(xspace)-1,1);
tao_0=0.07;
delta_a_hat=0.2;
a_hat=0.8;
coef=1/(tao_0*delta_a_hat);
for i=1:length(xspace)-1
 psi(i)=coef*func_g(t/tao_0,(xspace(i)-a_hat)/delta_a_hat);
end

