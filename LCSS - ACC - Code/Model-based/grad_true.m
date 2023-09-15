function grad = grad_true(A,B,Q,R,K)

nx=size(A,1);

%First we compute the term E_K

%Compute PK

PK=Q;
n_samples=5e3;
for i=1:n_samples
    PK = Q + K'*R*K + (A-B*K)'*PK*(A-B*K) ;
end

EK=(R+B'*PK*B)*K - B'*PK*A;

%Compute the covariance of the state Sigma_K

t=1000;
N_samples=1e4;
Sigma_s=[];
for j=1:N_samples
    x=normrnd(0,1,[nx,1]);
    Sigma=zeros(nx,nx);
    for i=1:t
        Sigma=Sigma + x*x';
        x=(A-B*K)*x;
    end
    Sigma_s{j}=Sigma;
end

%Averaging the covariance matrix

Sigma_K=zeros(nx,nx);

for i=1:N_samples
    Sigma_K = Sigma_K + (1/N_samples)*Sigma_s{i};
end

%True gradient 

grad=2*EK*Sigma_K;

end 