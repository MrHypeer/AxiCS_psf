% PoissonFct.m
% 
% Generation of a random number based on the Poisson distribution, with
% parameter lambda. 
% Based on the implementation proposed in "Numerical Recipes in C"
%
% Input:
% lambda        input noiseless value 
% 
% Output:
% S             output value with Poisson noise
%               same size as lambda input
%
% Copyright (C) 2015-2017  Nicolas Pavillon, Osaka University

function S = PoissonFct(lambda)

S = zeros(length(lambda),1);

for i = 1:length(lambda)
    if (lambda(i) < 12)
        %For small lambdas, direct method. We add random waiting times up
        %to reaching lambda.
        k = 0;
        produ = rand;
        lim = exp(-lambda(i));
        while produ >= lim
            produ = produ*rand;
            k = k+1;
        end
        S(i) = k;
    else
        %Rejection method for larger lambda. Poisson is bell-shaped, so
        %use a Lorentzian as comparison function.
        
        %Pre-computation of some values
        sq = sqrt(2*lambda(i));
        alxm = log(lambda(i));
        g = lambda(i)*alxm-gammaln(lambda(i)+1);
        
        firstT = 1;
        while (firstT)
            firstE = 1;
            while firstE %Find
                y = tan(pi*rand); %Random lorentzian value
                em = sq*y+lambda(i); %Scaling to ensure larger than our distribution
                if (em >= 0) %Security against negative tan
                    firstE = 0;
                end
            end
            em = floor(em); %Take integer as Poisson is discrete
            %Ratio to comparison function to find value (correction as our
            %random value is uniform)
            %Lorentzian is \frac{1}{\pi} \left( \frac{1}{1+y^2} \right)
            t = 0.9*(1+y*y)*exp(em*alxm-gammaln(em+1)-g); 
            
            if (t >= rand) %Rejection criterion
                firstT = 0;
            end
        end
        S(i) = em;
    end
end

end
