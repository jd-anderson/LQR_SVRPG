function grad = ZO2P(r,ns,A,B,Q,R,Kc,U)


nx=size(A,1);
nu=size(B,2);


K1=[];
K2=[];
cost_emp1=[];
cost_emp2=[];
%U=[];
%x=normrnd(0,0.25,[nx,1]);
for i=1:ns

    %U{i}=sample_matrix(nu,nx,r);
    
    %Sample policies
    K1{i}=Kc+U{i};
    K2{i}=Kc-U{i};

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

grad=zeros(nu,nx);
for i=1:ns
    grad=grad + ((nx*nu)/(2*ns*(r^2)))*(cost_emp1{i}-cost_emp2{i})*U{i};
end


end