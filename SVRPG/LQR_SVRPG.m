function [costs,grad_norm] = LQR_SVRPG(A,B,Q,R,K0,N,T,r,ns,ns_hat,eta,initial_gap,cost_optimal,cost_initial,q,U,U2)


costs=(cost_initial-cost_optimal)/initial_gap;
K=[];
K{1}=K0;
K_tilde=K{1};
grad_norm=[];
for n=1:N

    %Full gradient estimation
    g_full=ZO2P(1e-4,ns,A,B,Q,R,K_tilde,U);
    eta=eta_temp;
    for t=1:T
      
    
      %SVRG correction term

      ct=correction_term(r,ns_hat,A,B,Q,R,K_tilde,K{t},q,U2);
      
      %Stochastic gradient estimation

      v_t=g_full+ct;
      
      grad_norm=[grad_norm norm(v_t)];

      %Policy gradient update

      K{t+1}=K{t}-eta*(v_t);
      costs=[costs (cost(A,B,Q,R,K{t+1})-cost_optimal)/initial_gap];

    end
    
    Ktemp=K;
    K=[];
    K{1}=Ktemp{T+1};%initializing the next epoch
    K_tilde=K{1};%snapshot of the current policy   
end
costs=costs(:,1:length(costs)-1);
grad_norm=grad_norm(:,1:length(grad_norm)-1);
end