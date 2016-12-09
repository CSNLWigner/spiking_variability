function v = triuvals(a)
    v = a(find(~tril(ones(size(a)))));
end