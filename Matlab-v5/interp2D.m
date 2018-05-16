function [ yVal ] = interp2D( x1Vec, x2Vec, yMat, x1Val, x2Val )
% Function that does a dwo dimensionsional linear interpolation
% Matlab has its own in-built one of these - that has a different syntax 

%Check that both x1Vec and x2Vec are n_i * 1 vectors
sizex1Vec = size(x1Vec);
sizex2Vec = size(x2Vec);

if (sizex1Vec(1)<=1) || (sizex2Vec(1)<=1) || (sizex1Vec(2)~=1) || (sizex2Vec(2)~=1)
    error('Both x1Vec and x2Vec need to be n_i by 1 vectors. The length of the vectors may differ')
end

%Check whether x1Vec and x2Vec are monotonic
isx1VecMonotone = all(diff(x1Vec)>=0);
isx2VecMonotone = all(diff(x2Vec)>=0);

if (isx1VecMonotone*isx2VecMonotone == 0)
    error('One of x1Vec or x2Vec is non-monotonic');
end

n1 = sizex1Vec(1);

%Find where x1val fits into x1Vec and pick out the vectors in yMat
%conditional on these values
if (x1Val<=x1Vec(1))
    x1LowerIndex = 1; 
    x1UpperIndex = 2;        
elseif (x1Val>=x1Vec(n1))
    x1UpperIndex = n1;
    x1LowerIndex = n1 - 1;
else %x1Val is between the lower and the upper components of x1Vec
    %Get a vect
    x1VecLessX1Val = x1Vec - x1Val;
    x1VecLessX1Val(x1VecLessX1Val <0) = NaN;
    [nearestVal x1UpperIndex] = min(x1VecLessX1Val);
    x1LowerIndex = x1UpperIndex - 1;
end

    newX1VecShort = NaN(2, 1);
    newX2VecShort = NaN(2, 1);
    newYMatAfterInterpOnX2 = NaN(2, 1);
    
    yMatCondOnx1Lower = yMat(x1LowerIndex, :);
    yMatCondOnx1Upper = yMat(x1UpperIndex, :);
    
    newYMatAfterInterpOnX2(1) = interp1(x2Vec,yMatCondOnx1Lower,x2Val, 'linear', 'extrap');
    newYMatAfterInterpOnX2(2) = interp1(x2Vec,yMatCondOnx1Upper,x2Val, 'linear', 'extrap');
   
    newX1VecShort(1) = x1Vec(x1LowerIndex);
    newX1VecShort(2) = x1Vec(x1UpperIndex);
    
    yVal =  interp1(newX1VecShort,newYMatAfterInterpOnX2,x1Val, 'linear', 'extrap');
    
end

