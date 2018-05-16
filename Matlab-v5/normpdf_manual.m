function [ pdf] = normpdf_manual(x, mu, sigma)

pdf = ((sqrt(2.*pi*(sigma^2))).^-1) * (exp( - ((x - mu).^2)./(2*(sigma^2))));

end

