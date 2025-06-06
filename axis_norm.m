function norm_value = axis_norm(x,y,value,axis_lower,axis_upper,ax)
%This function takes in a value, such as the reduced velocity, determines
%where that value is in normalized position on the axis, then provides the
%normalized value of that position
ax_offset = get(ax,'Position');
ax_limits = get(ax,'XLim');
crossingIndices = find(y(1:end-1) < value & y(2:end) >= value);

% Initialize an array to store crossing points
crossingXValues = [];

% Interpolate to find the exact crossing point(s)
for i = 1:length(crossingIndices)
    idx = crossingIndices(i);
    % Linear interpolation
    x1 = x(idx);
    x2 = x(idx + 1);
    y1 = y(idx);
    y2 = y(idx + 1);
    
    % Calculate the x-value where the line crosses the threshold
    xCross = x1 + (value - y1) * (x2 - x1) / (y2 - y1);
    crossingXValues = [crossingXValues, xCross];
end

axis_upper = axis_upper-axis_lower;
xCross = xCross-axis_lower;
norm_value = ax_offset(1)+(xCross-axis_lower)/(axis_upper-axis_lower)*ax_offset(3);
end