# this module hosts helper functions to simplify notations in the big functions. 
# It's a try

# Production functions

function production_func(k,l,aalpha::Real,markup=1.0)
    y = k.^(aalpha).*l.^(1.0 .-aalpha)
    w = (1 .-aalpha).*markup.*(y./l)
    pi_per_k = (aalpha .* markup .+ (1.0 .-markup)) .* (y./k)
    return(y=y, w=w, pi_per_k=pi_per_k)
end

# here, use multiple dispatch
function production_func(k,l, aalpha,markup=1.0)
    output = get_y(k,l, aalpha;markup=markup)
    return(output)
end

# this function switches two variables, used in calc_equilibrium_and_update
function switch_two_var(a, b)
    temp = a
    a = b
    b=temp
    return(a,b)
end

# the Following returns utility function (as written by Moritz)
# return: util, util_c_deriv, labor_part 
# so far, the function only supports when inputs and outputs are scalars 
# For calling this function, input: c (number), labor (number), iii (index, integer)
# and last input is mpr_economy (our custom defined object)
# if you only want one output among the three, you can use util = util_func(..all inputs).util or .labor_part
# to extract the relevant fields

function util_fun(consumption, labor, iii, thtbar_vec, ies_vec, tht  ; sqrt_eps = sqrt(eps(Float64)))
    thtbar = thtbar_vec[iii]
    ies  = ies_vec[iii]

    if ((1.0 - sqrt_eps) < ies) & (ies < (1.0 + sqrt_eps))
        labor_part = 1.0
        util = log(consumption) - thtbar * tht/(1.0 + tht) * labor^((1.0+tht)/tht)
        util_c_deriv = 1.0/consumption
    else 
        labor_part = ( 1.0 + (1.0/ies - 1.0) * thtbar * tht/(1.0 + tht) * labor^((1.0+tht)/tht) )^(1.0/ies)
        util = (consumption^(1.0 - 1.0/ies)) * labor_part
        util_c_deriv = util/consumption
    end
    return (util = util, util_c_deriv = util_c_deriv, labor_part = labor_part)
end


function movsum(a, n::Int)
    out = similar(a, length(a))
    m = Int(floor(n/2)+1)
    out[1] = sum(a[1:m])
    for i =2:m
        out[i] = out[i-1]+a[i]
    end
    for i = m+1:length(a)-m+1
        out[i] = out[i-1]+a[i+m-1]-a[i-m]
    end
    for i =length(a)-m+2:length(a)
        out[i] = out[i-1]-a[i-m]
    end
    return out
end


function ols_detrend(y::Array{Float64,1})
    n = length(y)
    X = ones(n,2);
    X[:,2] = (1:n);
    return X*((X'*X)\X'*y)
end


function get_quadrature_points(n_grid::Int64,mu::Real=0.0, sig::Real=1.0) 
    nodes, weights = gausshermite(n_grid);
    nodes = nodes*sig*sqrt(2)+mu*ones(n_grid);
    weights = weights/sqrt(pi);
    return nodes, weights;
end


# the function below implements the Brent's method

# The following is the main function
# Inputs: x_a, x_b and a function f
# output: x_s, the root
# it has two different methods (when x_d is defined and when x_d is not defined yet)

function brent_solve(f::Function,x_a,x_b;max_iter=1000, tol=5E-16)
    y_a = f(x_a); y_b = f(x_b)
    
    if abs(y_a) < abs(y_b)
        
        x_a,x_b = switch_two_var(x_a,x_b)
        y_a,y_b = switch_two_var(y_a,y_b)
    end
    mflag = true
    x_c = x_a; y_c = y_a #initialize x_c
    x_s=Inf; y_s = Inf
    x_d = Inf #initialize x_d

    iter =0
    while (iter < max_iter) & (abs(y_s) > tol) & (abs(x_a - x_b) > tol)
        if y_a*y_b >0
            error("root not between a and b")
        end
        iter = iter +1 
        println("iteration number $iter")
       # println("current f(s) is $y_s")
       # println("x_a,x_b,y_a,y_b are $x_a, $x_b,$y_a,$y_b")
        # println("mflag is $mflag")
        
        x_s, mflag = brent_new_guess(x_a, x_b, x_c,y_a, y_b,y_c,mflag; x_d=x_d)
        # calc f(s)
        y_s = f(x_s)
        # assign new variables
        x_d=x_c
        x_c=x_b
        y_c=y_b
       # println("x_a,x_b are now $x_a, $x_b")
       # println("x_c, y_c are $x_c, $y_c")
        if (y_a*y_s < 0.0)
            x_b = x_s
            y_b = y_s
        else
            x_a = x_s
            y_a = y_s
        end
        
        if abs(y_a) < abs(y_b)
            x_a,x_b = switch_two_var(x_a,x_b)
            y_a,y_b = switch_two_var(y_a,y_b)
        end
        

    end
    return(x_s)

end





# the following gives a new guess x_s; can be used in other parts
function brent_new_guess(x_a,x_b,x_c,y_a,y_b,y_c,mflag;δ=1E-12,x_d=Inf)
    
    # obtain updated position
    if  (abs(x_b-x_c)>1E-6) & (abs(x_b-x_c)>1E-6)  # inverse quad interp
        x_s =  inv_quad_interp(x_a,x_b,x_c, y_a, y_b, y_c) 
    else # secant method
        x_s = x_b - y_b * (x_b - x_a)/(y_b - y_a)
       # println("x_s is $x_s")
    end   
    
    # println("x_s now is $x_s")

    # Next, check many conditions to decide whether to use the bisection method
    if isinf(x_d)
        cond = bisect_cond(x_a,x_b,x_c,x_s,mflag;δ=δ)
    else
        cond = bisect_cond(x_a,x_b,x_c,x_d,x_s,mflag;δ=δ)
    end

    if cond
        x_s = (x_a + x_b)/2
        mflag = true
    else
        mflag = false 
    end
    # print("mflag is $mflag")
    return(x_s = x_s, mflag = mflag)
end

# use multiple dispatch, get the sum
function inv_quad_interp(x_a,x_b,x_c,y_a,y_b,y_c)
    x_1 = inv_quad_interp(x_a, y_a, y_b, y_c)
    x_2 = inv_quad_interp(x_b, y_b, y_a, y_c)
    x_3 = inv_quad_interp(x_c, y_c, y_a, y_b)
    return(x_1 + x_2 + x_3)
    
end

# element wise interp
function inv_quad_interp(x_a, y_a, y_b, y_c)
    output = x_a * y_b * y_c / ((y_a-y_b)*(y_a-y_c))
    return(output)
end

# use multiple dispatch to solve the indicator for whether to use bisection
function bisect_cond(x_a,x_b,x_c,x_s,mflag; δ=1E-12)
    cond = false
    # first condition whether x_s is NOT between (3*x_a + x_b) and x_b
    cond1 = ((3*x_a + x_b)/4 < x_s < x_b) | (x_b < x_s < (3*x_a + x_b)/4)
    cond1 = !cond1
    cond = cond | cond1

    if mflag # mflag is a boolean variable
        cond2= abs(x_s-x_b) >= abs(x_b - x_c)/2
        cond3 = abs(x_b - x_c) < abs(δ)
        cond = cond | cond2 | cond3
    end
    
    return(cond)
end

# there are extra conditions to check if x_d is already defined
function bisect_cond(x_a,x_b,x_c,x_d,x_s,mflag;δ=1E-4)
    cond = bisect_cond(x_a,x_b,x_c,x_s,mflag; δ=δ) # first, grab the condition without x_d
    # print("cond = $cond")
    if !mflag # if mflag is false
        cond1 = abs(x_s - x_b) >= abs(x_b - x_d)/2
        # println("cond1 = $cond1")
        cond2 = abs(x_b - x_d) < 1E-12
        # println("cond2 = $cond2")

        cond = cond | cond1 | cond2
    end
    
    return(cond)    
end

function switch_two_var(a, b)
    temp = a
    a = b
    b=temp
    return(a,b)
end