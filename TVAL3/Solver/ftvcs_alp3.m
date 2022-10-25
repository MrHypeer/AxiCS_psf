% ftvcs_alp3.m
% 
% This function is a direct generalization of the original ftvcs_alp
% solver, for 3D matrices
%
% Copyright (C) 2015-2017  Nicolas Pavillon, Osaka University

function [U, out] = ftvcs_alp3(A,b,p,q,r,opts)
%
% Goal: solve   min sum ||D_i u||    (with or without the constraint u>=0) 
%                  s.t. Au = b
%       to recover image/signal u from encoded b,
%       which is equivalent to solve       min sum ||w_i||
%                                          s.t. D_i u = w_i
%                                               Au = b
% ftvcs_al solves Augmented Lagrangian function:
% 
% min_{u,w} sum ||w_i|| - sigma'(Du-w) - delta'(Au-b) 
%                   + beta/2 ||Du-w||_2^2 + mu/2||Au-b||_2^2 ,
%
% by an alternating algorithm:
% i)  while norm(up-u)/norm(up) > tol_inn
%     1) Fix w^k, do Gradient Descent to 
%            - sigma'(Du-w^k) - delta'(Au-b) + beta/2||Du-w^k||^2 + mu/2||Au-f||^2;
%            u^k+1 is determined in the following way:
%         a) compute step length tau > 0 by BB formula
%         b) determine u^k+1 by
%                  u^k+1 = u^k - alpha*g^k,
%            where g^k = -D'sigma - A'delta + beta D'(Du^k - w^k) + mu A'(Au^k-f), 
%            and alpha is determined by Amijo-like nonmonotone line search;
%     2) Given u^k+1, compute w^k+1 by shrinkage
%                 w^k+1 = shrink(Du^k+1-sigma/beta, 1/beta);
%     end
% ii) update Lagrangian multipliers by
%             sigma^k+1 = sigma^k - beta(Du^k+1 - w^k+1)
%             delta^k+1 = delta^k - mu(Au^k+1 - b).
% iii)accept current u as the initial guess to run the loop again
%
% Inputs:
%       A        : either an matrix representing the measurement or a struct 
%                  with 2 function handles:
%                           A(x,1) defines @(x) A*x;
%                           A(x,2) defines @(x) A'*x;
%       b        :  either real or complex input vector representing the
%                   noisy observation of a grayscale image
%       p, q     :  size of original image
%       opts     :  structure to restore parameters
%
%
% variables in this code:
%
% lam1 = sum ||wi||
% lam2 = ||Du-w||^2 (at current w).
% lam3 = ||Au-f||^2
% lam4 = sigma'(Du-w)
% lam5 = delta'(Au-b)
%
%   f  = lam1 + beta/2 lam2 + mu/2 lam3 - lam4 - lam5
%
%   g  = A'(Au-f)
%   g2 = D'(Du-w) (coefficients beta and mu are not included)
%
%
% Numerical tests illustrates that this solver doestn't require large beta
% and mu. ( <100 usually)
%
%
% Written by: Chengbo Li
% Advisor: Prof. Yin Zhang and Wotao Yin
% Computational and Applied Mathematics department, Rice University
% May. 2, 2009

%% Initializations and setting options

global D Dt
[D,Dt] = defDDt3;

% problem dimension
n = p*q*r;

% unify implementation of A
if ~isa(A,'function_handle')
    A = @(u,mode) f_handleA(A,u,mode); 
end

% get or check opts
opts = ftvcs_al_opts(opts); 

% mark important constants
mu = opts.mu;
beta = opts.beta;
tol_inn = opts.tol_inn;
tol_out = opts.tol;
gam = opts.gam;

%% Scaling

% check if A*A'=I
tmp = rand(length(b),1);
if norm(A(A(tmp,2),1)-tmp,1)/norm(tmp,1) < 1e-3
    opts.scale_A = false;
end
clear tmp;

% check scaling A (scaling is done by eigenvalues)
if opts.scale_A
    [mu,A,b] = ScaleA(n,mu,A,b,opts.consist_mu);
end 

% check scaling b (scales b between 0.5-1.5, and scales mu accordingly)
if opts.scale_b
    [mu,b,scl] = Scaleb(mu,b,opts.consist_mu);
end

% calculate A'*b
Atb = A(b,2);

%% Prepare working variables

% initialize U, beta
muf = mu;
betaf = beta;     % final beta
[U,mu,beta] = ftvcs_al_init(p,q,r,Atb,scl,opts);    % U: p*q
if mu > muf; mu = muf; end
if beta > betaf; beta = betaf; end
muDbeta = mu/beta;
rcdU = U;

% initialize multiplers
sigmax = zeros(p,q,r);                       % sigmax, sigmay: p*q 
sigmay = zeros(p,q,r);
sigmaz = zeros(p,q,r);
delta = zeros(length(b),1);                % delta: m

% initialize D^T sigma + A^T delta
DtsAtd = zeros(p*q*r,1); 

% initialize out.n2re
if isfield(opts,'Ut')
    Ut = opts.Ut*scl;        %true U, just for computing the error
    nrmUt = normfro(Ut);
else
    Ut = []; 
end
if ~isempty(Ut)
    out.n2re = normfro(U - Ut)/nrmUt; 
end

% prepare for iterations
out.mus = mu; out.betas = beta;
out.res = []; out.itrs = []; out.f = []; out.obj = []; out.reer = [];
out.lam1 = []; out.lam2 = []; out.lam3 = []; out.lam4 = []; out.lam5 = [];
out.itr = Inf;
out.tau = []; out.alpha = []; out.C = []; gp = [];
out.cnt = [];

%% First computation of Lagrangian function

%D an Dt are defined in defDDt, which computes differential matrices
[Ux,Uy,Uz] = D(U);                   % Ux, Uy: p*q
if opts.TVnorm == 1
    Wx = max(abs(Ux) - 1/beta, 0).*sign(Ux);
    Wy = max(abs(Uy) - 1/beta, 0).*sign(Uy);
    Wz = max(abs(Uz) - 1/beta, 0).*sign(Uz);
    lam1 = sum(sum(sum(abs(Wx) + abs(Wy) + abs(Wz))));
else
    V = sqrt(Ux.*conj(Ux) + Uy.*conj(Uy) + Uz.*conj(Uz));        % V: p*q
    V(V==0) = 1;
    S = max(V - 1/beta, 0)./V;        % S: p*q
    Wx = S.*Ux;                       % Wx, Wy: p*q
    Wy = S.*Uy;
    Wz = S.*Uz;
    lam1 = sum(sum(sum(sqrt(Wx.*conj(Wx) + Wy.*conj(Wy) + Wz.*conj(Wz)))));  
end  

%Get intermediate values for Lagrangian function
[lam2,lam3,lam4,lam5,f,g2,Au,g] = get_g(U,Ux,Uy,Uz,Wx,Wy,Wz,lam1,beta,mu,A,b,...
    Atb,sigmax,sigmay,sigmaz,delta);
%lam, f: constant      g2: pq        Au: m         g: pq

% compute gradient
d = g2 + muDbeta*g - DtsAtd;

%% Start outer loop: gradient descent
count = 1;
Q = 1; C = f;                     % Q, C: costant
out.f = [out.f; f]; out.C = [out.C; C];
out.lam1 = [out.lam1; lam1]; out.lam2 = [out.lam2; lam2]; out.lam3 = [out.lam3; lam3];
out.lam4 = [out.lam4; lam4]; out.lam5 = [out.lam5; lam5];

for ii = 1:opts.maxit
    if opts.disp
        fprintf('outer iter = %d, total iter = %d, normU = %4.2e; \n',count,ii,normfro(U));
    end
    
    % compute tau first
    if ~isempty(gp)
        dg = g - gp;                        % dg: pq
        dg2 = g2 - g2p;                     % dg2: pq
        ss = uup'*uup;                      % ss: constant
        sy = uup'*(dg2 + muDbeta*dg);       % sy: constant
        % sy = uup'*((dg2 + g2) + muDbeta*(dg + g));
        % compute BB step length
        tau = abs(ss/max(sy,eps));               % tau: constant
        
        clear uup; %Not used during smaller loops 
        fst_itr = false;
    else
        % do Steepest Descent at the 1st ieration
        %d = g2 + muDbeta*g - DtsAtd;         % d: pq
        [dx,dy,dz] = D(reshape(d,p,q,r));                    %dx, dy: p*q
        dDd = normfro(dx)^2 + normfro(dy)^2 + normfro(dz)^2;      % dDd: cosntant
        clear dx dy dz;
        Ad = A(d,1);                        %Ad: m
        % compute Steepest Descent step length
        tau = abs((d'*d)/(dDd + muDbeta* (Ad'*Ad) ));
        clear Ad;
        
        % mark the first iteration 
        fst_itr = true;
    end    
    
    % keep the previous values
    Up = U; gp = g; g2p = g2; Aup = Au; Uxp = Ux; Uyp = Uy; Uzp = Uz;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ONE-STEP GRADIENT DESCENT %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    taud = tau*d;
    U = U(:) - taud;
    % projected gradient method for nonnegtivity
    if opts.nonneg
        U = max(real(U),0);
    end
    U = reshape(U,p,q,r);                    % U: p*q (still)
    [Ux,Uy,Uz] = D(U);                        % Ux, Uy: p*q
    
    [lam2,lam3,lam4,lam5,f,g2,Au,g] = get_g(U,Ux,Uy,Uz,Wx,Wy,Wz,lam1,beta,mu,A,b,...
        Atb,sigmax,sigmay,sigmaz,delta);
   
    % Nonmonotone Line Search
    alpha = 1;
    du = U - Up;                          % du: p*q
    const = opts.c*beta*(d'*taud);
	clear taud;
	
   %% Start inner loop: steepest descent
    % Unew = Up + alpha*(U - Up)
    cnt = 0; flag = true; maxcnt = opts.maxcnt;
    while f > C - alpha*const
        if cnt == maxcnt
            % shrink gam
            gam = opts.rate_gam*gam;

            % give up and take Steepest Descent step
            if opts.disp
                disp(['    count of back tracking attains ', int2str(maxcnt)]);
            end

            %d = g2p + muDbeta*gp - DtsAtd;
            [dx,dy,dz] = D(reshape(d,p,q,r));
            dDd = normfro(dx)^2 + normfro(dy)^2 + normfro(dz)^2;
            clear dx dy dz;
            Ad = A(d,1);
            tau = abs((d'*d)/(dDd + muDbeta* (Ad'*Ad)));
            clear Ad;
            U = Up(:) - tau*d;
            % projected gradient method for nonnegtivity
            if opts.nonneg
                U = max(real(U),0);
            end
            U = reshape(U,p,q,r);
            [Ux,Uy,Uz] = D(U);
            if opts.TVnorm == 1
                % ONE-DIMENSIONAL SHRINKAGE STEP
                aux = Ux - sigmax/beta;
                Wx = max(abs(aux) - 1/beta, 0).*sign(aux);
                aux = Uy - sigmay/beta;
                Wy = max(abs(aux) - 1/beta, 0).*sign(aux);
                aux = Uz - sigmaz/beta;
                Wz = max(abs(aux) - 1/beta, 0).*sign(aux);
                lam1 = sum(sum(sum(abs(Wx) + abs(Wy) + abs(Wz))));
                clear aux;
            else
                % TWO-DIMENSIONAL SHRINKAGE STEP
                Uxbar = Ux - sigmax/beta;
                Uybar = Uy - sigmay/beta;
                Uzbar = Uz - sigmaz/beta;
                V = sqrt(Uxbar.*conj(Uxbar) + Uybar.*conj(Uybar) + Uzbar.*conj(Uzbar)); % V: p*q
                V(V==0) = 1;
                S = max(V - 1/beta, 0)./V;                         % S: p*q
                clear V;
                Wx = S.*Uxbar;
                Wy = S.*Uybar;
                Wz = S.*Uzbar;
                clear Uxbar Uybar Uzbar;
                lam1 = sum(sum(sum(sqrt(Wx.*conj(Wx) + Wy.*conj(Wy) + Wz.*conj(Wz)))));
            end
            
            [lam2,lam3,lam4,lam5,f,g2,Au,g] = get_g(U,Ux,Uy,Uz,Wx,Wy,Wz,lam1,...
                beta,mu,A,b,Atb,sigmax,sigmay,sigmaz,delta);
            alpha = 0; % remark the failure of back tracking
            break;
        end
        if flag
            dg = g - gp;
            dg2 = g2 - g2p;
            dAu = Au - Aup;                 % dAu: m
            dUx = Ux - Uxp;
            dUy = Uy - Uyp;
            dUz = Uz - Uzp;
            flag = false;
        end
        alpha = alpha*opts.gamma;
        [lam3,lam5,Au,g,g2] = update_g1(alpha,gp,dg,g2p,dg2,Aup,dAu,b,delta);
        U = Up + alpha*reshape(du,p,q,r);
        [lam2,lam4,f,Ux,Uy,Uz] = update_g2(lam1,lam3,lam5,alpha,beta,mu,...
                Wx,Wy,Wz,Uxp,dUx,Uyp,dUy,Uzp,dUz,sigmax,sigmay,sigmaz);
        cnt = cnt + 1;
    end
    
    % if back tracking is succeceful, then recompute
    if alpha ~= 0
        if opts.TVnorm == 1
            % ONE-DIMENSIONAL SHRINKAGE STEP
            aux = Ux - sigmax/beta;
            Wx = max(abs(aux) - 1/beta, 0).*sign(aux);
            aux = Uy - sigmay/beta;
            Wy = max(abs(aux) - 1/beta, 0).*sign(aux);
            aux = Uz - sigmaz/beta;
            Wz = max(abs(aux) - 1/beta, 0).*sign(aux);
            clear aux;
        else
            % TWO-DIMENSIONAL SHRINKAGE STEP
            Uxbar = Ux - sigmax/beta;
            Uybar = Uy - sigmay/beta;
            Uzbar = Uz - sigmaz/beta;
            V = sqrt(Uxbar.*conj(Uxbar) + Uybar.*conj(Uybar) + Uzbar.*conj(Uzbar));
            V(V==0) = 1;
            S = max(V - 1/beta, 0)./V;
            clear V;
            Wx = S.*Uxbar;
            Wy = S.*Uybar;
            Wz = S.*Uzbar;
            clear Uxbar Uybar Uzbar;
        end
        
        
        % update parameters related to Wx, Wy
        [lam1,lam2,lam4,f,g2] = update_W(beta,...
            Wx,Wy,Wz,Ux,Uy,Uz,sigmax,sigmay,sigmaz,lam1,lam2,lam4,f,opts.TVnorm);
    end
    
    % update reference value
    Qp = Q; Q = gam*Qp + 1; C = (gam*Qp*C + f)/Q;
    uup = U - Up; uup = uup(:);           % uup: pq
    nrmuup = norm(uup,'fro');                   % nrmuup: constant
    
    out.res = [out.res; nrmuup];
    out.f = [out.f; f]; out.C = [out.C; C]; out.cnt = [out.cnt;cnt];
    out.lam1 = [out.lam1; lam1]; out.lam2 = [out.lam2; lam2]; out.lam3 = [out.lam3; lam3];
    out.lam4 = [out.lam4; lam4]; out.lam5 = [out.lam5; lam5];
    out.tau = [out.tau; tau]; out.alpha = [out.alpha; alpha];

    if ~isempty(Ut), out.n2re = [out.n2re; norm(U - Ut,'fro')/norm(Ut,'fro')]; end

    % compute gradient
    d = g2 + muDbeta*g - DtsAtd;

    nrmup = normfro(Up);
    RelChg = nrmuup/nrmup;

    %% Stopping conditions
    if RelChg < tol_inn && ~fst_itr
        if (opts.disp)
            disp('tol_inn satisfied');
        end
        count = count + 1;
        RelChgOut = normfro(U-rcdU)/nrmup;
        out.reer = [out.reer; RelChgOut];
        rcdU = U;
        out.obj = [out.obj; f + lam4 + lam5];
        if isempty(out.itrs)
            out.itrs = ii;
        else
            out.itrs = [out.itrs; ii - sum(out.itrs)];
        end

        % stop if already reached final multipliers
        if RelChgOut < tol_out || count > opts.maxcnt 
            if (opts.disp)
                disp('tol_out satisfied');
            end
            if opts.isreal
                U = real(U);
            end
            if exist('scl','var')
                U = U/scl;
            end
            out.itr = ii;
            fprintf('Number of total iterations is %d. \n',out.itr);
            return
        else
            if (opts.disp)
                fprintf('   RelChgOut %.2e\n', RelChgOut);
            end
        end
        
        % update multipliers
        [sigmax,sigmay,sigmaz,delta,lam4,lam5,~] = update_mlp(beta,mu, ...
            Wx,Wy,Wz,Ux,Uy,Uz,Au,b,sigmax,sigmay,sigmaz,delta,lam4,lam5,f);
 
        % update penality parameters for continuation scheme
        beta0 = beta;
        beta = beta*opts.rate_ctn;
        mu = mu*opts.rate_ctn;
        if beta > betaf; beta = betaf; end
        if mu > muf; mu = muf; end
        muDbeta = mu/beta;
        out.mus = [out.mus; mu]; out.betas = [out.betas; beta];

        % update function value, gradient, and relavent constant
        f = lam1 + beta/2*lam2 + mu/2*lam3 - lam4 - lam5;
        DtsAtd = -(beta0/beta)*d;     % DtsAtd should be divded by new beta instead of the old one for consistency!  
        d = g2 + muDbeta*g - DtsAtd;
        
        %initialize the constants
        gp = [];
        gam = opts.gam; Q = 1; C = f;
    else
        if (opts.disp)
            fprintf('   RelChg %.2e\n', RelChg);
        end
    end

end

%% Final adjustments for ouptut
if opts.isreal
    U = real(U);
end
if exist('scl','var')
    fprintf('Attain the maximum of iterations %d. \n',opts.maxit);
    U = U/scl;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%END OF ftvcs_alp3

function [lam2,lam3,lam4,lam5,f,g2,Au,g] = get_g(U,Ux,Uy,Uz,Wx,Wy,Wz,lam1,...
    beta,mu,A,b,Atb,sigmax,sigmay,sigmaz,delta)
global Dt

% A*u 
Au = A(U(:),1);

% g
g = A(Au,2) - Atb;

% lam3
Aub = Au-b;
lam3 = norm(Aub,'fro')^2;

%lam5
lam5 = delta'*Aub;
clear Aub;

% lam2
Vx = Ux - Wx;
Vy = Uy - Wy;
Vz = Uz - Wz;
lam2 = sum(sum(sum(Vx.*conj(Vx) + Vy.*conj(Vy) + Vz.*conj(Vz))));

% g2 = D'(Du-w)
g2 = Dt(Vx,Vy,Vz);

%lam4
lam4 = sum(sum(sum(conj(sigmax).*Vx + conj(sigmay).*Vy + conj(sigmaz).*Vz)));

% f
f = lam1 + beta/2*lam2 + mu/2*lam3 - lam4 - lam5;

function [lam3,lam5,Au,g,g2] = update_g1(alpha,gp,dg,g2p,dg2,Aup,dAu,b,delta)
Au = Aup + alpha*dAu;
Aub = Au-b;
lam3 = norm(Aub,'fro')^2;
lam5 = delta'*Aub;

g = gp + alpha*dg;
g2 = g2p + alpha*dg2;


function [lam2,lam4,f,Ux,Uy,Uz] = update_g2(lam1,lam3,lam5,...
    alpha,beta,mu,Wx,Wy,Wz,Uxp,dUx,Uyp,dUy,Uzp,dUz,sigmax,sigmay,sigmaz)

Ux = Uxp + alpha*dUx;
Uy = Uyp + alpha*dUy;
Uz = Uzp + alpha*dUz;

aux = Ux - Wx;
lam2 = sum(sum(sum(aux.*conj(aux))));
lam4 = sum(sum(sum(conj(sigmax).*aux)));
aux = Uy - Wy;
lam2 = lam2 + sum(sum(sum(aux.*conj(aux))));
lam4 = lam4 + sum(sum(sum(conj(sigmay).*aux)));
aux = Uz - Wz;
lam2 = lam2 + sum(sum(sum(aux.*conj(aux))));
lam4 = lam4 + sum(sum(sum(conj(sigmaz).*aux)));

f = lam1 + beta/2*lam2 + mu/2*lam3 - lam4 - lam5;



function [lam1,lam2,lam4,f,g2] = update_W(beta,...
    Wx,Wy,Wz,Ux,Uy,Uz,sigmax,sigmay,sigmaz,lam1,lam2,lam4,f,option)
global Dt

% update parameters because Wx, Wy were updated
tmpf = f -lam1 - beta/2*lam2 + lam4;
if option == 1
    lam1 = sum(sum(sum(abs(Wx) + abs(Wy) + abs(Wz))));
else
    lam1 = sum(sum(sum(sqrt(Wx.^2 + Wy.^2 + Wz.^2))));
end
Vx = Ux - Wx;
Vy = Uy - Wy;
Vz = Uz - Wz;
g2 = Dt(Vx,Vy,Vz);
lam2 = sum(sum(sum(Vx.*conj(Vx) + Vy.*conj(Vy) + Vz.*conj(Vz))));
lam4 = sum(sum(sum(conj(sigmax).*Vx + conj(sigmay).*Vy + conj(sigmaz).*Vz)));
f = tmpf +lam1 + beta/2*lam2 - lam4;



function [sigmax,sigmay,sigmaz,delta,lam4,lam5,f] = update_mlp(beta,mu, ...
    Wx,Wy,Wz,Ux,Uy,Uz,Au,b,sigmax,sigmay,sigmaz,delta,lam4,lam5,f)

tmpf = f + lam4 + lam5;

Aub = Au-b;
lam5 = delta'*Aub;
delta = delta - mu*Aub;
clear Aub;

aux = Ux - Wx;
sigmax = sigmax - beta*aux;
lam4 = sum(sum(sum(conj(sigmax).*aux)));
aux = Uy - Wy;
sigmay = sigmay - beta*aux;
lam4 = lam4 + sum(sum(sum(conj(sigmay).*aux)));
aux = Uz - Wz;
sigmaz = sigmaz - beta*aux;
lam4 = lam4 + sum(sum(sum(conj(sigmaz).*aux)));

f = tmpf - lam4 - lam5;


function [U,mu,beta] = ftvcs_al_init(p,q,r,Atb,scl,opts)

% initialize mu beta
if isfield(opts,'mu0')
    mu = opts.mu0;
else
    error('Initial mu is not provided.');
end
if isfield(opts,'beta0')
    beta = opts.beta0;
else
    error('Initial beta is not provided.');
end

% initialize U
%Opts.init is a field for initial guess
[mm,nn,oo] = size(opts.init); 
if max( [mm,nn,oo] ) == 1
	%If no guess provided, either use zero, or inverted observations as first guess
    switch opts.init
        case 0, U = zeros(p,q,r);
        case 1, U = reshape(Atb,p,q,r);
    end
else
    U = opts.init*scl;  
    if mm ~= p || nn ~= q || oo ~= r
        fprintf('Input initial guess has incompatible size! Switch to the default initial guess. \n');
        U = reshape(Atb,p,q,r);
    end
end

