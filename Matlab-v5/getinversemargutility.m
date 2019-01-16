function [ invmargut ] = getinversemargutility (margut)

% Whats going on here?

global gamma

if gamma == 1
    invmargut = 1./margut;
else
    invmargut = margut.^(-1/gamma);
end

end

