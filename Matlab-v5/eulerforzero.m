function [euler  ] = eulerforzero(A0, A1, Y)

%%
%-------------------------------------------------------------------------------%
%This function returns the following quantity:
%u'(c_t) - b(1+r)u'(c_t+1) 
%This quantity =0 when the Euler equation u'(c_t) = b(1+r)u'(c_t+1)  
%is satified with equality
%-------------------------------------------------------------------------------%

%% ------------------------------------------------------------------------ 
% Declare global we need this file have access to

global r beta interpMethod 
global Edu1 linEdU1 linearise
global Agrid1

%% ------------------------------------------------------------------------ 
% Get marginal utility 
%Get marginal utility at consumption tomorrow. 
if linearise == 0
    du1AtA1 = interp1(Agrid1,Edu1,A1, interpMethod, 'extrap');
elseif linearise == 1
    invDu1atA1 = interp1(Agrid1,linEdU1,A1, interpMethod, 'extrap');
    du1AtA1 = getmargutility(invDu1atA1);
end    
    
%% ------------------------------------------------------------------------ 
% Check whether tomorrow's (expected) marginal utility negative If so throw an error

if (du1AtA1 < 0)    
   error('approximated marginal utility in negative')
end
      
%% ------------------------------------------------------------------------ 
% Get consumption today and the required output
todaycons = A0 + Y - A1/(1+r);
euler = getmargutility(todaycons) - (beta * (1+r) * du1AtA1) ;

end

