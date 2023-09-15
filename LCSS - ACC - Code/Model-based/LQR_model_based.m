function [costs,grad_norm] = LQR_model_based(A,B,Q,R,K_0,N,eta,initial_gap,cost_optimal,cost_initial)


costs=(cost_initial-cost_optimal)/initial_gap;
K=[];
K{1}=K_0;
grad_norm=[];

for n=1:N-1

    grad=grad_true(A,B,Q,R,K{n});
    grad_norm=[grad_norm norm(grad,"fro")];
    K{n+1}=K{n}-eta*grad;
    costs=[costs (cost(A,B,Q,R,K{n+1})-cost_optimal)/initial_gap];
    
end


end 

