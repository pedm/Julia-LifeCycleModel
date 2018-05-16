function [ margut ] = getmargutility (cons)

global gamma

if cons<=0
   error('Consumption is <=0');
end
                                
if gamma == 1
    margut = 1./(cons);
else
    margut = (cons).^(-gamma);
end

end

