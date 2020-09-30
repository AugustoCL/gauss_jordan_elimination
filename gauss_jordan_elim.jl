using DelimitedFiles, LinearAlgebra

function troca_linha(i::Int64, nlinha::Int64)
    for n = (i+1):nlinha        # iterate over lines above to check if could be swap
        if A[n,i] != 0          # condition to swap row
            L = copy(A[i,:])    # copy line to swap
            A[i,:] = A[n,:]     # swap occur
            A[n,:] = L
            break
        end
    end
end

function gauss_jordan(A::Array)::Array
    # check if matrix is singular
    LinearAlgebra.det(A[:,1:end-1]) == 0 ? throw("Insert a non-singular matrix") : nothing

    nlinha = size(A)[1]
    for i = 1:nlinha
        if A[i,i] == 0.0                            # check if need rows swap
            troca_linha(i, nlinha)
        end

        A[i,:] = A[i,:]./ A[i,i]                    # divide pivot line by pivot element

        for j = 1:nlinha                            # iterate each line for each pivot column, except pivot line
            if j != i                               # jump pivot line
                A[j,:] = A[j,:] - A[i,:]*A[j,i]     # apply gauss jordan in each line
            end
        end
    end
    return A
end

A = readdlm("input.txt")
gauss_jordan(A)
A = readdlm("input_singular_matrix.txt")
gauss_jordan(A)

A = [1 2 3; 4 5 6]
gauss_jordan(A)
A = [1 2 3; 4 8 12]
gauss_jordan(A)

### logic in the for loops
# A[1,:] = A[1,:]./ A[1,1]
# A[2,:] = A[2,:] - A[1,:]*A[2,1]
# A[3,:] = A[3,:] - A[1,:]*A[3,1]

# A[2,:] = A[2,:]./ A[2,2]
# A[1,:] = A[1,:] - A[2,:].*A[1,2]
# A[3,:] = A[3,:] - A[2,:].*A[3,2]

# A[3,:] = A[3,:]./ A[3,3]
# A[1,:] = A[1,:] - A[3,:].*A[1,3]
# A[2,:] = A[2,:] - A[3,:].*A[2,3]