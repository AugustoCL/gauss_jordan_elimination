using LinearAlgebra, BenchmarkTools


# Gauss-Jordan elimination functions ------------------------------------------------

function swap_rows(i::T, nlinha::T) where {T<:Integer}
    for n ∈ (i+1):nlinha        # iterate over lines above to check if could be swap
        if A[n,i] ≠ 0.0         # condition to swap row
            L = copy(A[i,:])    # copy line to swap
            A[i,:] = A[n,:]     # swap occur
            A[n,:] = L
            break
        end
    end
end

function gauss_jordan(A::Matrix{T}) where {T<:Number}
    
    # convert to float to avoid InexactError: Int64()
    (T <: Integer) && (A = convert.(Float64, A))

    # check if matrix is singular
    m, n = size(A)
    if m == n
        @assert det(A) ≠ 0.0 "Must insert a non-singular matrix"
    else
        @assert det(A[:,1:end-1]) ≠ 0.0 "Must insert a non-singular matrix or a system matrix [A b]"
    end

    for i ∈ axes(A, 1)
        if A[i,i] == 0.0                            # check if need swap rows
            swap_rows(i, m)
        end

        @. A[i,:] = A[i,:] / A[i,i]                 # divide pivot line by pivot element

        for j ∈ axes(A, 1)                          # iterate each line for each pivot column, except pivot line
            if j ≠ i                                # jump pivot line
                @. A[j,:] = A[j,:] - A[i,:]*A[j,i]  # apply gauss jordan in each line
            end
        end
    end

    return A
end


# Applying function in some matrix ----------------------------------------------------

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
@btime gauss_jordan2(A)

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


## logic in the for loops -------------------------------------------------
# A[1,:] = A[1,:]./ A[1,1]
# A[2,:] = A[2,:] - A[1,:]*A[2,1]
# A[3,:] = A[3,:] - A[1,:]*A[3,1]

# A[2,:] = A[2,:]./ A[2,2]
# A[1,:] = A[1,:] - A[2,:].*A[1,2]
# A[3,:] = A[3,:] - A[2,:].*A[3,2]

# A[3,:] = A[3,:]./ A[3,3]
# A[1,:] = A[1,:] - A[3,:].*A[1,3]
# A[2,:] = A[2,:] - A[3,:].*A[2,3]