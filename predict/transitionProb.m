function A = transitionProb(X , n)

% Return an estimate of the Transition Probabilities matrix given the
% Markov Chain X.
%
% sum(A , 2) = 1 , A = Pr(x_{k} | x_{k - 1} );
%
% M = 100;
% K = 1000;
% X = ceil(M*rand(1 , K));
% A = Proba_A(X);

if ( nargin < 2 | isempty(n) )
   
    T = full(sparse(X(1:end - 1) , X(2 : end) , 1));
   
else
   
    T = full(sparse(X(1:end - 1) , X(2 : end) , 1) , n , n);
   
end

sum_T = sum(T , 1);

A = (T./sum_T(ones(size(T , 1) , 1) , :))';