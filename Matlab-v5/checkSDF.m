function checkSDF(policyA1, EV, dU)

global T beta r
global interpMethod Agrid numPointsA numPointsY


%Create a matrix that stores whether one is borrowing constrained or not
borrowingConstrained = NaN(T,numPointsA, numPointsY);
sdf = NaN(T, numPointsA, numPointsY);

%The below could be re-written in matrix notation for a substantial speed-up

h = 0.000000001;
for ixt=T:-1:1 %Loop from last time period to the first one
  tomAssGrid = Agrid(ixt + 1, :);              
  lbA1 = tomAssGrid(1);         
  for ixA = 1:1:numPointsA %Loop through each point on the asset grid                             
        for ixY = 1:1:numPointsY      
            
            A1 = policyA1(ixt, ixA, ixY);
            
            %Are we Borrowing Constrained
                if (A1 < lbA1 + 0.0001)
                    borrowingConstrained(ixt, ixA, ixY) = 1;
                else
                    borrowingConstrained(ixt, ixA, ixY) = 0;                    
                end
            
            tomExValFunc = EV(ixt + 1, :, ixY);
            upper = interp1(tomAssGrid,tomExValFunc, A1 + h,interpMethod,'extrap');
            lower = interp1(tomAssGrid,tomExValFunc, A1 - h,interpMethod,'extrap');
            deriv = (upper - lower)/(2*h);
            sdf(ixt, ixA, ixY) =  dU(ixt, ixA, ixY) / ((1+ r) * beta * deriv);            
        end
    end
end

%When we are borrowing constrained put a -1 in the stochastic discount
%factor matrix
    sdf(borrowingConstrained == 1) = -1;

%Now count the sdfs that are greater that     
    %How many non borrowing constrained points in the state space are there
    nonBC = T*numPointsA*numPointsY - sum(sum(sum(borrowingConstrained)));
    
    tolSDF = NaN(4, 1);
    tolSDF(1, 1) = 0.25;
    tolSDF(2, 1) = 0.1;
    tolSDF(3, 1) = 0.05;
    tolSDF(4, 1) = 0.025;
    
    %How many SDFs are outside these tolerances
    for ixTol = 1:1:4
        failSDF = NaN(T, numPointsA, numPointsY);
        lowerTol = (1 - tolSDF(ixTol, 1));
        upperTol = (1 + tolSDF(ixTol, 1));
        failSDF( sdf < lowerTol | sdf > upperTol) = 1;
        failSDF( sdf > lowerTol & sdf < upperTol) = 0;
        failSDF( borrowingConstrained == 1) = 0;
        failRep = sum(sum(sum(failSDF))) / nonBC;
        fprintf('Prop. of SDFs outside range %5.3f and %5.3f is: %5.3f.\n',lowerTol, upperTol, failRep)
    end
     
    
end
