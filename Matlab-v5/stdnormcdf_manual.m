function [ cdf] = stdnormcdf_manual(x)

%This function gives the cdf of a standard normal

approxforminusinf = - 20;
cdf = quad(@(t) stdnormpdf_manual(t), approxforminusinf, x, 1e-12);

end

