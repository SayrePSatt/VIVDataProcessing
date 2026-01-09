function [X,Y] = uniqueXYExtraction(x,y)
if (y(2)-y(1)) == 0
    kk = 1;
    for i = 1:length(y)
        X(kk,1) = x(i);
        if (y(i+1)-y(i)) ~= 0
            break
        end
        kk = kk+1;
    end

    Y(1,1) = y(1);
    kk = 2;
    for i = 2:length(y)
        if (y(i)-y(i-1)) ~= 0
            Y(kk,1) = y(i);
            kk = kk+1;
        end
    end
else
    kk = 1;
    for i = 1:length(x)
        Y(kk,1) = y(i);
        if (x(i+1)-x(i)) ~= 0
            break
        end
        kk = kk+1;
    end

    X(1,1) = x(1);
    kk = 2;
    for i = 2:length(x)
        if (x(i)-x(i-1)) ~= 0
            X(kk,1) = x(i);
            kk = kk+1;
        end
    end
end
end