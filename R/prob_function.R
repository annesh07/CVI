f0 = function(x){tail(cumsum(x),1)-cumsum(x)}
f1 = function(x){cumsum(x)*(tail(cumsum(x),1)-cumsum(x))}
