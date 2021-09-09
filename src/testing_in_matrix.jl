using BenchmarkTools

include("src/gauss_jordan_elim.jl");

# regular matrix
A = [1 2 3
     4 5 6];
gauss_jordan(A)

# a singular matrix
A = [1 2 3 
     4 8 12];
gauss_jordan(A)

# less performatic because convert to float to avoid InexactError: Int64()
A = [1 4 8  1  4
     4 5 6  8  11
     1 3 2  4  8
     4 5 67 23 0];
@btime gauss_jordan(A)

# more performatic because it's already Float64
A = Float64[1 4 8  1  4
            4 5 6  8  11
            1 3 2  4  8
            4 5 67 23 0];
@btime gauss_jordan(A)

A = Float64[1 4 8  1
            4 5 6  8
            1 3 2  4
            4 5 67 23];
gauss_jordan(A)

# when dimensions doesn't match
A = [1 4 8 1  4;
     4 5 6 8 11;
     1 3 2 4  8];
gauss_jordan(A)