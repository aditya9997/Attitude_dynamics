k = 4;
S = zeros(k, k);
for n = 1:k
    for m = 0:(n - 1)
        if m == 0
            if n == 1
                S(n,m + 1) = 1;
            else
                S(n,m + 1) = S(n - 1, m + 1) * ( (2 * n - 1 ) / n );
            end
        elseif m == 1
            S(n, m + 1) = S(n, m) * sqrt( (n - m + 1) * 2 / (n + m) );
        else
            S(n, m + 1) = S(n, m) * sqrt( (n - m + 1) / (n + m) );
        end
    end
end

S = S';
S