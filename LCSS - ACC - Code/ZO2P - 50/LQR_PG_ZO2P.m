function [costs,grad_norm] = LQR_PG_ZO2P(A,B,Q,R,K_0,N,r,ns,eta,initial_gap,cost_optimal,cost_initial,U)


costs=(cost_initial-cost_optimal)/initial_gap;
K=[];
K{1}=K_0;
grad_norm=[];
for n=1:N-1

    %Zeroth-order gradient estimation with one-point estimation

    grad=ZO2P(r,ns,A,B,Q,R,K{n},U);
    grad_norm=[grad_norm norm(grad)];
    
    %Policy gradient update

    K{n+1}=K{n}-eta*grad;
    
    costs=[costs (cost(A,B,Q,R,K{n+1})-cost_optimal)/initial_gap];

end


end 

