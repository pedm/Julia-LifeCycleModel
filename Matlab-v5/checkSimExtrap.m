function [ a1] = checkSimExtrap( lba1,y )

% This is a function to check, in selecting next period's asset we haven't
% selected a value that is less than the borrowing constraint. This could
% occur for one of two reasons. First it could be that we are extrapolating
% (i.e. income in the period is larger than the largest income in the grid.
% If this is the case we set next period's assets to the lowest permissable
% level. Otherwise it is likely that there is an error - and in this case
% we cause the programme to stop

global Ygrid numPointsY

if (y > Ygrid(numPointsY)) || (y < Ygrid(1)) 
    a1 = lba1;
else
    error('Next periods asset is below minimum permissable assets. And we are not extrapolating. Likely there is a bug')
end


end

