function [ x] = stdnorminv_manual(p)

%This routines gives the standard normal value associated with a particular p value,
%with truncations at -3 and 3

if (p<stdnormcdf_manual(-3))
    x = -3;
elseif  (p>stdnormcdf_manual(3))
    x = 3;
else
    boundforzero = [-3, 3];
    x = fzero(@(t) stdnormcdf_manual(t) - p, boundforzero, optimset('TolX',1e-12));
end 

end

