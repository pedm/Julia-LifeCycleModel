function [ pdf] = stdnormpdf_manual(x)

%This function gives the pdf of a standard normal

pdf = ((sqrt(2.*pi)).^-1) * (exp( - ((x).^2)/(2)));

end

