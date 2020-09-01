# Example wth julia Package SigmoidalProgramming
# Install Sigmoidal Programming by entering pkg mode (type ']'  in the repl)
# then type 'add https://github.com/madeleineudell/SigmoidalProgramming.jl'
using SigmoidalProgramming
using Random
using DataStructures

function SP()
    # generate problem data
    Random.seed!(5)
    nvar = 200
    nineqconstr = 20
    l = -rand(nvar)
    u = rand(nvar)
    A = rand(nineqconstr, nvar)
    b = rand(nineqconstr)
    z = zeros(nvar)
    fs = fill(logistic, nvar)
    dfs = fill(logistic_prime, nvar)
    problem = LinearSP(fs, dfs, z, A, b)

    # branch and bound to solve the problem
    # pq is a priority queue of the branch and bound nodes at the leaves of the tree
    # bestnodes is a list of the best branch and bound nodes found, in the order they were found
    pq, bestnodes, lbs, ubs = solve_sp(l, u, problem)

    # the best node found yet is the top node on the priority queue
    node = dequeue!(pq)
    # println("best node has node.ub = $(node.ub) and solution $(node.x)")

    # lbs and ubs record the upper and lower bounds on the optimal value
    # found at each iteration
    println("lbs: ",lbs)
    println("ubs: ",ubs)
    nothing
end

@time SP()
@code_warntype SP()
