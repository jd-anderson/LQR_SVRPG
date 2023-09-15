function random_matrix = sample_matrix(m,n,r)
    
    random_points = normrnd(0,1,[m, n]);
    random_matrix = r*(random_points/norm(random_points,"fro"));
    
end