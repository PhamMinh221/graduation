max = 0;
for i = 1 : 10000
    if max < norm((Q((i+1)*T)-Q(i*T))*d)
        max = norm((Q((i+1)*T)-Q(i*T))*d);
    end
end
