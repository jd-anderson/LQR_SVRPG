function cost = cost_oracle(A,B,K,R,Q,x)

cost=0;
tau=1e4;
for t=1:tau
    cost=cost + (x'*Q*x + (K*x)'*R*(K*x));
    x=(A-B*K)*x;
end

end