# MRP IP: an integer program that ensures the production of each SKU occurs as late as possible to minimize inventory costs
# Note: This will be ultimately be modeled under Dantzig-Wolfe Decomposition

# ---- Abstract Formulation ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

# ---- Concrete Formulation ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

P := 3;
T := 2;
K := 2;

minimize Earliness: x[1,1] + x[2,1] + x[3,1] + x[1,2] + x[2,2] + x[3,2]                                                                       ;

s.t. Demand11     : x[1,1] + x[2,1] + x[3,1]                                                                                            >= b[1] ;
s.t. Demand21     : x[1,1] + x[2,1] + x[3,1]                                                                                            >= b[2] ;
s.t. Demand31     : x[1,1] + x[2,1] + x[3,1]                                                                                            >= b[3] ;
s.t. Demand12     : x[1,1] + x[2,1] + x[3,1] + x[1,2] + x[2,2] + x[3,2]                                                               >= b[4] ;
s.t. Demand22     : x[1,1] + x[2,1] + x[3,1] + x[1,2] + x[2,2] + x[3,2]                                                               >= b[5] ;
s.t. Demand32     : x[1,1] + x[2,1] + x[3,1] + x[1,2] + x[2,2] + x[3,2]                                                               >= b[6] ;

s.t. Resource11   : x[1,1] + x[2,1] + x[3,1]                                                                                            <= b[7] ;
s.t. Resource21   : x[1,1] + x[2,1] + x[3,1]                                                                                            <= b[8] ;
s.t. Resource12   : x[1,1] + x[2,1] + x[3,1] + x[1,2] + x[2,2] + x[3,2]                                                               <= b[9] ;
s.t. Resource22   : x[1,1] + x[2,1] + x[3,1] + x[1,2] + x[2,2] + x[3,2]                                                               <= b[10];

s.t. LotSize11    : x[1,1]                                                   + y[1,1]                                                    >= b[11];
s.t. LotSize21    :           x[2,1]                                                   + y[2,1]                                          >= b[12];
s.t. LotSize31    :                     x[3,1]                                                   + y[3,1]                                >= b[13];
s.t. LotSize12    :                               x[1,2]                                                   + y[1,2]                      >= b[14];
s.t. LotSize22    :                                         x[2,2]                                                   + y[2,2]            >= b[15];
s.t. LotSize32    :                                                   x[3,2]                                                   + y[3,2]  >= b[16];

s.t. Production11 : x[1,1]                                                   + y[1,1]                                                    >= b[17];
s.t. Production21 :           x[2,1]                                                   + y[2,1]                                          >= b[18];
s.t. Production31 :                     x[3,1]                                                   + y[3,1]                                >= b[19];
s.t. Production12 :                               x[1,2]                                                   + y[1,2]                      >= b[20];
s.t. Production22 :                                         x[2,2]                                                   + y[2,2]            >= b[21];
s.t. Production32 :                                                   x[3,2]                                                   + y[3,2]  >= b[22];

# ---- DW Formulation ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# ---- Decomposition Description ----

# 0. add slack variables to all constraints to change all inequalities to equalities

# 1. there will be a master problem to satisfy the 'LotSize' and 'Production' constraints

# 2. there will be a subproblem for each time period to satisfy the 'Demand' and 'Resource' constriants that apply to those time periods

    # in the following order:
    # subproblem 1 will satisfy time period 1 constraints during time period 0
	# subproblem 2 will satisfy time period 1 & 2 constraints during time period 0, using the solution found in subproblem 1 as an initial basis
	# ...
	# subproblem T will satisfy time period 1, 2, ... , T constraints during time period 0, using the solution found in subproblem (T - 1) as an initial basis
	
	# when the solution from a subproblem is fed into the successor subproblem as an initial basis, primal feasibility must be evaluated and potentially solved via dual simplex

# 3. when the optimal solution is found, it will be used to set production requirements for time period 1

# 4. time period 1 will then become time period 0, time period 2 becomes time period 1, ... , time period T becomes time period (T - 1)

    # when time period 1 becomes time period 0, its production requirements must correspond to the updates made in the data file (ie. less resources are now available)
	
# 5. the solution process resets and uses this optimal solution as a starting basis to find the next optimal solution

# ---- Decomposition Structure ----

# 1 Master Problem
# 2 subproblems

# ---- Master Problem: LotSize & Production constraints ----

set J1;    # set of extreme points of Subproblem 1
set J2;    # set of extreme points of Subproblem 2
set K1;    # set of extreme rays of Subproblem 1
set K2;    # set of extreme rays of Subproblem 2

param x1{J1};    # the extreme points representing the amount of SKU i produced at time t for Subproblem 1
param x2{J2};    # the extreme points representing the amount of SKU i produced at time t for Subproblem 2
param w1{K1};    # the extreme rays representing the amount of SKU i produced at time t for Subproblem 1
param w2{K2};    # the extreme rays representing the amount of SKU i produced at time t for Subproblem 2
param c1;    # the objective function vector of coefficients (ie. Earliness) for Subproblem 1
param c2;    # the objective function vector of coefficients (ie. Earliness) for Subproblem 2
param D1;    # the left hand side matrix of coefficients representing LotSize and Production constraints for Subproblem 1
param D2;    # the left hand side matrix of coefficients representing LotSize and Production constraints for Subproblem 2
param b0;    # the right hand side vector of values representing LotSize and Production constraints

var lam1{J1} >= 0;    # the vector of values to create a convex combination of extreme points for Subproblem 1
var lam2{J2} >= 0;    # the vector of values to create a convex combination of extreme points for Subproblem 2
var thea1{K1} >= 0;    # the vector of values to create a linear combination of extreme rays for Subproblem 1
var thea2{K2} >= 0;    # the vector of values to create a linear combination of extreme rays for Subproblem 2

minimize Z0: sum{j in J1}(lam1[j] * c1 * x1[j]) + sum{k in K1}(thea1[k] * c1 * w1[k]) + sum{j in J2}(lam2[j] * c2 * x2[j]) + sum{k in K2}(thea2[k] * c2 * w2[k]);    # ensure that production of each SKU occurs as late as possible
s.t. Master: sum{j in J1}(lam1[j] * D1 * x1[j]) + sum{k in K1}(thea1[k] * D1 * w1[k]) + sum{j in J2}(lam2[j] * D2 * x2[j]) + sum{k in K2}(thea2[k] * D2 * w2[k]) = b0;    # satisy all LotSize and Production constraints: read 'Abstract Formulation' above for each of those descriptions
s.t. Convex1: sum{j in J1}(lam1[j]) = 1;    # the sum decison variables representing a set of convex values, for Subproblem 1, must equal 1: this is what makes it a convex combination as opposed to a linear combination
s.t. Convex2: sum{j in J2}(lam2[j]) = 1;    # the sum decison variables representing a set of convex values, for Subproblem 2, must equal 1: this is what makes it a convex combination as opposed to a linear combination

# ---- Subproblem 1: Demand & Resource constraints across time period 1 ----

param q;    # the dual variable vector for constraint 'Master' in 'Master Problem'
param r1;    # the dual variable vector for constraint 'Convex1' in 'Master Problem' 
param c1;    # the objective function vector of coefficients (ie. Earliness) for time period 1
param D1;    # the left hand side matrix of coefficients representing LotSize and Production constraints for time period 1
param F1;    # the left hand side matrix of coefficients representing Demand and Resource constraints for time period 1
param b1;    # the right hand side vector of values representing Demand and Resource constraints for time period 1

var X1 >= 0;    # the extreme point in the polyhedron for Subproblem 1 

minimize Z1: ((c1 - (q * D1)) * X1) - r1;    # minimize the reduced cost of the extreme point in the polyhedron for Subproblem 1
s.t. Sub1: (F1 * X1) = b1;    # satisy Demand and Resource constraints for time period 1: read 'Abstract Formulation' above for each of those descriptions

# ---- Subproblem 2: Demand & Resource constraints across time periods 1 & 2 ----

param q;    # the dual variable vector for constraint 'Master' in 'Master Problem'
param r2;    # the dual variable vector for constraint 'Convex2' in 'Master Problem' 
param c2;    # the objective function vector of coefficients (ie. Earliness) for time periods 1 & 2
param D2;    # the left hand side matrix of coefficients representing LotSize and Production constraints for time periods 1 & 2
param F2;    # the left hand side matrix of coefficients representing Demand and Resource constraints for time periods 1 & 2
param b2;    # the right hand side vector of values representing Demand and Resource constraints for time periods 1 & 2

var X2 >= 0;    # the extreme point in the polyhedron for Subproblem 2 # initialize with basis found in Subproblem 1

minimize Z2: ((c2 - (q * D2)) * X2) - r2;    # minimize the reduced cost of the extreme point in the polyhedron for Subproblem 2
s.t. Sub2: (F2 * X2) = b2;    # satisy Demand and Resource constraints for time periods 1 & 2: read 'Abstract Formulation' above for each of those descriptions
