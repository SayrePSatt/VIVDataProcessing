function str = num2roman(x)
    % Validates input is a scalar integer 
    assert(isscalar(x) && floor(x)==x, 'Input must be a scalar integer');
    assert(1 <= x && x <= 3999, 'Input must be between 1 and 3999');
    
    numbers = [1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1];
    letters = {'m', 'cm', 'd', 'cd', 'c', 'xc', 'l', 'xl', 'x', 'ix', 'v', 'iv', 'i'};
    str = '';
    num = x;
    
    for i=1:numel(numbers)
        while (num >= numbers(i))
            str = [str letters{i}];
            num = num - numbers(i);
        end
    end
end