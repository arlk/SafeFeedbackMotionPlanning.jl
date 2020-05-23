function _chebpts(N)
    K = N-1
    n = (cos.(π*(K:-1:0)/K) .+ 1)./2
    w = similar(n)

    Kh = K % 2 == 0 ? K^2 - 1 : K^2
    w[1] = w[end] = 0.5/Kh

    Kt = div(K, 2)
    for k = 1:Kt
        wk = 0.0
        for j = 0:Kt
            β = (j == 0) || (j == K/2) ? 1 : 2
            wk += β/K/(1-4*j^2)*cos(2*π*j*k/K)
        end
        w[K-k+1] = w[k+1] = wk
    end

    return n, w
end

function _chebpoly(nodes, D)
    # D : Max degree

    N = length(nodes)
    T = zeros(D, N)

    for i = 1:D
        for j = 1:N
            if i == 1
                T[i,j] = 1
            elseif i == 2
                T[i,j] = nodes[j]
            else
                T[i,j] = 2*nodes[j]*T[i-1,j] - T[i-2,j]
            end
        end
    end

    return T
end

function _chebpolyder(T, nodes, D)
    # D : Max degree

    N = length(nodes)
    dT = similar(T)

    for i = 1:D
        for j = 1:N
            if i == 1
                dT[i,j] = 0
            elseif i == 2
                dT[i,j] = 1
            else
                dT[i,j] =  2*T[i-1,j] + 2*nodes[j]*dT[i-1,j] - dT[i-2,j]
            end
        end
    end

    return dT
end
