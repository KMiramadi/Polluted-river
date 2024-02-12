%%Create matrix A
%0<y<1
clc
N=1000;
h_y=1/N;
eta=0.01; %
beta=0.2;
tao_0=0.07; %ALso in psimaker
a_hat=0.8; %Also in psimaker
a_roof=sqrt(1/4*func_w(0)^2+beta^2)-1/(2*func_w(0));
delta_a_hat=0.2; %Also in psimaker
yspace=0:h_y:1;
T_end=2;


A=zeros(length(yspace)-1,length(yspace)-1);
A(1,2)=func_a_j(eta,h_y,yspace(1))+func_c_j(eta,h_y,yspace(1));
A(1,1)=func_b_j(eta,h_y,yspace(1))-(2*h_y*a_roof*func_a_j(eta,h_y,yspace(1))/eta);

for row=2:N-1
      A(row,row-1)=func_a_j(eta,h_y,yspace(row));
      A(row,row)=func_b_j(eta,h_y,yspace(row));
      A(row,row+1)=func_c_j(eta,h_y,yspace(row));
end

A(end,end)=func_b_j(eta,h_y,yspace(N));
A(end,end-1)=func_a_j(eta,h_y,yspace(N));



%%
h_t=1/1000;
T_end=1.5;
tspace=0:h_t:T_end;

%U=zeros(length(yspace),length(tspace));

time_init_cond=zeros(length(yspace)-1,1);

[t, sol]=CrankNicolNew(A,h_y,h_t,T_end,time_init_cond);
%%
dirchlet_end=zeros(1,length(tspace));
sol=cat(1,sol,dirchlet_end);


%%
figure
mesh(t,yspace,sol)
xlabel('Time')
ylabel('y')
shading interp

%%
w=arrayfun(@(y) func_w(y) ,yspace);
plot(yspace,w,'LineWidth',1.5,'Color','b')
xlabel('y')
ylabel('width of river')
title('Plot of w(y)')
grid on
%%
u_div_w=zeros(length(yspace),length(t));
for i=1:length(t)
    for j=1:length(yspace)
        u_div_w(j,i)=sol(j,i)/w(j);
    end
   
end
%%
%Plot over time of u/w
figure
plot(yspace,u_div_w(:,180),'LineWidth',1.5)
xlabel('y')
title('u/w at \tau = 0.1790')
grid on
figure
plot(t,u_div_w(500,:),'color','g','LineWidth',1.5)
xlabel('\tau')
title('u/w as function of time at y=0.4990')
grid on

  

%%
%Part 3
%FIsh DEATH
%%
%u_infinity
gamma_inf=0.1;
gamma_tot=0.8;  
u_inf = @(y) gamma_inf.*max(u_div_w(round(y./h_y)+1,:))

full_u_inf=arrayfun( @(y) u_inf(y),yspace);
figure
plot(yspace,full_u_inf, 'LineWidth',1.5)
title('u_{infinity} as a function of y')
xlabel('y')
grid on
%%
%u_total
gamma_inf=0.1;
gamma_tot=0.8;  
u = @(tao,y) sol(round(y./h_y)+1,round(tao./h_t)+1);

u_tot = @(y) gamma_tot.*integral( @(tao) u(tao,y)./func_w(y),0,t(end));

full_total= arrayfun( @(y) u_tot(y), yspace);
figure
plot(yspace,full_total,'LineWidth',1.5)
title('u_{tot} as function of y')
grid on
%%
%sigma function
sigma = @(c) c./(1+c);
figure
plot(yspace,sigma(full_u_inf),'LineWidth', 1.5)
xlabel('y')
title('p_{\infty} as function of y')
grid on

figure
plot(yspace,sigma(full_total),'LineWidth', 1.5)
xlabel('y')
title('p_{tot} as function of y')
grid on
%%
% Plot all
figure
plot(yspace,full_u_inf, 'LineWidth',1.5)
hold on
plot(yspace,full_total,'LineWidth',1.5)
hold on
plot(yspace,sigma(full_u_inf),'LineWidth', 1.5)
hold on
plot(yspace,sigma(full_total),'LineWidth', 1.5)
title('Plot of u_{\infty}, u_{tot}, p_{\infty}, p_{tot}')
legend('u_{\infty}','u_{tot}','p_{\infty}','p_{tot}')
xlabel('y')
grid on
%%
%Total share of dead fish
%%
%Here AbsTol and RelTol are neccesary. The computations otherwise takes an
%excessive amount of time with Arrayvalued as True. Probably due to the
%integral in u_tot
tic
shareTot_numerator=integral(@(x) sigma(u_tot(x))*func_w(x),0,1,'Arrayvalued',1,'AbsTol',1e-5,'RelTol',1e-3);
toc
denom=integral( @(y) func_w(y),0,1);
result_total_share=shareTot_numerator/denom
%%
shareINF_numerator=integral(@(x) sigma(u_inf(x))*func_w(x),0,1,'Arrayvalued',1);

result_inf_share=shareINF_numerator/denom
