% Scaling and accuracy tests with numerical integrations done in parallel

clear; close all
addpath("functions","plotting")

plot_figs = 1;
save_data = 0;

% Test 1: error in K as a function of M. We use the 1-layer modon problem
% with (U,a,R,beta) = (1,1,inf,0) where we have an analytical result of
% K = j_11 where j_11 is the first non-zero root of the Bessel function
% J_1(x). Note: the error decreases exponentially until it reaches a value
% set by the error in the numerical integration of the Bessel function
% integral (see JJ_int.m).

% set parameters:

mu = 0;
lambda = 0;
M_max = 12;
tol = 1e-10;

% derived parameters:

j_11 = fzero(@(x) besselj(1,x),[3.8,3.9]);
M = 2:M_max;
K_M = 0*M;

% calculate K for a given M value:

for m = M

    disp([' - Finding K for M = ' num2str(m)])

    A = zeros(m);
    B = zeros(m);
    c = [1/4; zeros(m-1,1)];
    d = (-1).^(0:m-1)';
    
    for j = 0:m-1
        for k = 0:m-1
            A(j+1,k+1) = JJ_int(@(x) A_func(x,lambda,mu),j,k,tol);
            B(j+1,k+1) = JJ_int(@(x) B_func(x,lambda,mu),j,k,tol);
        end
    end
    
    [K2,~] = EVP(A,B,[mu*c c],d,1);
    [~,i] = min(abs(K2-j_11^2));
    K_M(m+1-M(1)) = sqrt(K2(i));

end

% plot error |K-j_11| on a logarithmic plot:

if plot_figs == 1
    figure;
    semilogy(M,abs(K_M-j_11)); grid; xlabel('M'); ylabel('|K-j_{11}|')
end

% Test 2: time taken as a function of number of layer, N. We use a N-layer
% modon problem with random mu and lambda in each layer and seek active
% layer solutions. Note: the time taken increases approximately linearly
% with the number of layers. The limiting step is the evaluation of the
% numerical integrals, in particular the diagonal entries. Therefore the
% time scales with the number of diagonal entries, T ~ O(NM).
%
% Note: This text runs in parallel by parallelising the calculation of A_ij
% and B_ij across up to C = M^2 cores.

% set parameters:

M = 8;          % number of term in expansion
N_max = 100;     % maximum number of layers
S = 5;          % number of repeats for each N, final T is an average
C = 32;         % number of cores/parallel processes

% derived parameters:

mu = rand(N_max,S);
lambda = rand(N_max,S);
N = 1:N_max;
K_N = zeros(N_max,N_max,S);
a_N = zeros(N_max,N_max,M,S);
T_N = zeros(N_max,2);

% calculate K for a given N value:

parpool(C);

for n = N
    
    for is = 1:S
    
        disp([' - Finding K for N = ' num2str(n) ', Run ' num2str(is) ' of ' num2str(S)])
        
        tic

        % calculate A and B matrices:

        A = zeros(n,n,M^2);
        B = zeros(n,n,M^2);

        parfor i = 0:M^2-1      % parallelised calculation of A and B terms
            j = mod(i,M); k = floor(i/M);
            A(:,:,i+1) = JJ_int(@(x) A_func(x,lambda(1:n,is),mu(1:n,is)),j,k);
            B(:,:,i+1) = JJ_int(@(x) B_func(x,lambda(1:n,is),mu(1:n,is)),j,k);
        end

        A = reshape(permute(reshape(A,[n n M M]),[1 3 2 4]),[n*M n*M]);
        B = reshape(permute(reshape(B,[n n M M]),[1 3 2 4]),[n*M n*M]);

        % calculate B_i, c_i and d_j:
        
        Bcell = cell(n,1);
        c = zeros(n*M,n+1);
        d = zeros(n*M,n);

        for i = 1:n
            Bcell{i} = kron(eye(M),diag(1:n == i))*B;
            c(:,1) = c(:,1) + mu(i)*kron(eye(M),diag(1:n == i))*[ones(n,1); zeros((M-1)*n,1)]/4;
            c(:,i+1) = kron(eye(M),diag(1:n == i))*[ones(n,1); zeros((M-1)*n,1)]/4;
            d(:,i) = kron((-1).^(0:M-1)',(1:n == i)');
        end

        T_N(n,1) = T_N(n,1) + toc;
        
        % solve eignevalue problem with initial guess K^2 = K_0^2:
        
        tic

        if n == 1   % guess K^2 and a_n based on results from previous n
            K_02 = 20*ones(n,1);
            a_0 = zeros(n*M,1);
        else
            K_02 = [K_N(1:n-1,n-1,is); K_N(n-1,n-1,is)].^2;
            a_0 = reshape([reshape(a_N(1:n-1,n-1,1:M,is),[n-1 M]); reshape(a_N(n-1,n-1,1:M,is),[1 M]);],[n*M 1]);
        end

        [Kn,an] = EVP_optim(A,Bcell,c,d,K_02,a_0);
        K_N(1:n,n,is) = sqrt(Kn);
        a_N(1:n,n,1:M,is) = reshape(an,[n M]);

        T_N(n,2) = T_N(n,2) + toc;
    
    end

end

T_N = T_N/S;    % average T_N over S runs

% plot runtime against N:

if plot_figs == 1
    figure;
    plot(N,T_N); grid; xlabel('N'); ylabel('Runtime (s)'); legend('Numerical Integrations','Eigenvalue Root Finding')
end

% save all data:

if save_data == 1
    save('test_data.mat')
end
