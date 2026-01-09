function [yPoints,uSlice] = extractSliceVelocities(XX,YY,Y,ux,xLocation)
% Extract x,y,u,v
for i = 1:length(xLocation)
    yPoints = Y;
    for j = 1:length(yPoints)
        xPoints(j) = xLocation(i);
    end
    uSlice(:,i) = interp2(XX,YY,ux,xPoints',yPoints);
end
end