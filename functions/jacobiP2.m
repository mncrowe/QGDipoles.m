function P = jacobiP2(n,a,b,x)
% Evaluates the Jacobi polynomial using a faster method that the Matlab function 'JacobiP'

    P = 0;
    for  s = 0:n
        P = P + 1/(factorial(s)*factorial(n+a-s)*factorial(b+s)*factorial(n-s))*((x-1)/2).^(n-s).*((x+1)/2).^s;
    end
    P = P*factorial(n+a)*factorial(n+b);
end