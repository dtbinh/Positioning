function p = p_k(x, k)

x_0 = x(1:end-k);
x_k = x((1+k):end);
mx_0 = mean(x_0);
mx_k = mean(x_k);

nominator = sum( (x_0-mx_0) .* (x_k-mx_k) );
denominator =  sqrt(sum((x_0-mx_0).^2)) * sqrt(sum((x_k-mx_k).^2));

p = nominator/denominator;

end