for mu = 3:13
    for beta =  3:13
        for tol = [1e-2 1e-3 1e-4 1e-5 1e-6]
            opts.mu = 2^mu;
            opts.beta = 2^beta;
            opts.tol = tol;
        end
    end
end