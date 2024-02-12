function [t,U] = CrankNicolNew(A,h_y, delta_t,T_end,init_cond)
%Crank Nicolson method for CE3HT23
%A the matrix and b the vector b(t) for du/dt = Au+ b(t)
% init_cond is a vector with inital condition. It is assumed that we solve
% from 0 to T_end
steps=T_end/delta_t;
U=zeros(length(init_cond),steps);
U(:,1)=init_cond;
t=0:delta_t:T_end;



for k=1:length(t)-1
    b_n=psimaker(t(k),h_y);
 
    b_np1=psimaker(t(k+1),h_y);

    mat=speye(length(init_cond),length(init_cond))-1/2*delta_t*A;
    solFor=(speye(length(init_cond),length(init_cond))+1/2*delta_t*A)*U(:,k)+1/2*delta_t*(b_n+b_np1);
    sol= mat \ solFor;
    U(:,k+1)=sol;
end


    
end

