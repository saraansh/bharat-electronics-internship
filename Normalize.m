function [X_normalized] = func_norm(X)

minx = min(X);

maxx = max(X);

X_normalized = (2 * (X - minx) / (maxx - minx)) - 1;