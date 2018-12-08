# MRP IP: an integer program that ensures the production of each SKU occurs as late as possible to minimize inventory costs

set P;    # number of SKU's
param T;    # total number of time periods
param K;    # total number of available resources

param l{P};    # lead time for each SKU
param r{P, P};    # internal demand: the number of SKU i required to make one unit of SKU j
param d{P, 1..T};    # external demand for SKU i in time period t
param v{P};    # initial inventory of SKU i at time period 0
param q{P};    # minimum lot size for each SKU
param u{1..K, 1..T};    # fraction of resource k available during time period t

param M;    # a large number

var x{P, 1..T};    # amount of SKU i produced at time period t
var y{P, 1..T} binary;    # production of SKU i does/doesn't occur at time period t

minimize Earliness: sum{i in P, t in 1..T}((T - t) * x[i,t]);    # ensure that production of each SKU occurs as late as possible
s.t. Demand{i in P, t in 1..T}: sum{k in 1..(t - l[i])}(x[i,k]) + v[i] - sum{k in 1..t}(d[i,k] + sum{j in P}(r[i,j] * x[j,k])) >= 0;    # the supply availble, up to each time period, must satisfy the total internal and external demand of each SKU
s.t. LotSize{i in P, t in 1..T}: x[i,t] >= y[i,t] * q[i];    # every production batch must satisfy minimum lot size requirements of each SKU
s.t. Production{i in P, t in 1..T}: y[i,t] >= (x[i,t] / M);    # the amount of SKU i produced must be zero if production for that SKU doesn't occur
s.t. Resource{k in 1..K, t in 1..T}: sum{i in P}(u[i,k] * x[i,t]) <= 1;    # the availbility of each resource in each period cannot be exceeded
