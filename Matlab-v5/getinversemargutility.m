function [ invmargut ] = getinversemargutility (margut)

global gamma

if gamma == 1
    invmargut = 1./margut;
else
    invmargut = margut.^(-1/gamma);
end

end

