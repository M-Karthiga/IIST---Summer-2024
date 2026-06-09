function [P, R] = computePR(re)
    if re == 1
        error('re = 1 is singular. Choose re ≠ 1.');
    end
    if re > 1
        beta = log(re + sqrt(re^2 - 1)) / (re * sqrt(re^2 - 1));
    else
        beta = acos(re) / (re * sqrt(1 - re^2));
    end

    P = (re^2 - 1) / (2 * (re^2 + 1));
    R = (3 * re^4 * (beta * (2 * re^2 - 1) - 1)) / ...
        (16 * pi * (re^2 - 1) * (re^2 + 1));
end


