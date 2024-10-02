function y = Fu(x, alpha)
    w = 1 / (1 - alpha);
    y = 1 + x - (1 + (x) ^ w) ^ (1 / w);
end