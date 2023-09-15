function corr_term = correction_term(r,ns,A,B,Q,R,K_tilde,Kc,q,U)


nx=size(A,1);
nu=size(B,2);



corr_term=zeros(nu,nx);
for j=1:q
    K1=[];
    K2=[];
    cost_emp1=[];
    cost_emp2=[];
    for i=1:ns
     
        %Sample policies
    
        K1{i}=Kc+U{i};
        K2{i}=K_tilde+U{i};
    
        if sum(abs(eig(A-B*K1{i}))>=1)>0
              stop=1;
        end
        
        if sum(abs(eig(A-B*K2{i}))>=1)>0
              stop=1;
        end
    
    
        %Initial state distrubiton 
    
        x=normrnd(0,1,[nx,1]);
    
        %Compute empirical cost
    
        cost_emp1{i}=cost_oracle(A,B,K1{i},R,Q,x);
        cost_emp2{i}=cost_oracle(A,B,K2{i},R,Q,x);
    end
    
    ct=zeros(nu,nx);
    for i=1:ns
        ct=ct + ((nx*nu)/(ns*(r^2)))*(cost_emp1{i}-cost_emp2{i})*U{i};
    end
    corr_term=corr_term+(1/q)*ct;
end

end