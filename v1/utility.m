function [ utils ] = utility (cons)

%This function takes consumption as an argument and returns utility. The
%utility functin is CRRA except that we add a very small number (eps) to
%consumption so that the computer can deal wi

global gamma 

if cons<=0 
   error('Error in utility. Consumption is <=0');
end
                                

if gamma == 1
    utils = log(cons);
else
    utils = ((cons)^(1-gamma)  )/(1-gamma);
end

end

