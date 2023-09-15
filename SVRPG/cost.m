function cost = cost(A,B,Q,R,K)

nx=size(A,1);


%Compute PK

PK=Q;
n_samples=2e3;
for i=1:n_samples
    PK = Q + K'*R*K + (A-B*K)'*PK*(A-B*K) ;
end


%Cost for a fixed initial condition

x_0=ones(nx,1);

cost=x_0'*PK*x_0;


end



