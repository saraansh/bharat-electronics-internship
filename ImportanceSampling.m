% Generate samples
N = 2000;
t = 1/N;
X_samples = rand(N,1);

% Probability distribution function 'p'
p_func = @(x) tan(degtorad(x));

% Probability distribution function 'q'
q_func = @(x) sin(degtorad(x));

% Evaluate for each sample
target = p_func(X_samples);
proposal = q_func(X_samples);

% Calculate importance weights p(x)/q(x)
weights = target./proposal;

% Resample, with replacement, according to importance weights
approximates = randsample(X_samples, N, true, weights);

% Plotting curve for p(x)
figure
subplot(3,2,1)
x = linspace(0, 1, N);
plot(x, p_func(x))
title('P(x)')
axis normal

% Plotting random samples X
subplot(3,2,2)
hist(X_samples, N)
title('X or F(x)')
axis normal

% Plotting curve for q(x)
subplot(3,2,3)
plot(x, q_func(x))
title('Q(x)')
axis normal

% Plotting approximated X
subplot(3,2,4)
hist(approximates, N)
title('Approximated X')
axis normal

% Plotting curve for target values
subplot(3, 2, 5)
plot(x, sort(target))
title('Expected Curve')
axis normal

% Plotting curve for approximated values
subplot(3, 2, 6)
plot(x, p_func(sort(approximates)))
title('Approximated Curve')
axis normal

% Print sum for comparison
s1 = sum(target)
s2 = sum(p_func(approximates))

% Mean absolute error
error = mae(sort(target), p_func(sort(approximates)))