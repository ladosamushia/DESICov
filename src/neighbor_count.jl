using NearestNeighbors
using Base.Threads

@inline minimg(δ::Float32, L::Float32) = δ - round(δ / L) * L
@inline P2(μ::Float32) = muladd(1.5f0*μ, μ, -0.5f0)  # (3/2)μ^2 - 1/2
@inline function binidx(edges::Vector{Float32}, x::Float32)
    k = searchsortedlast(edges, x)
    k == 0 ? 0 : (k == length(edges) ? length(edges)-1 : k)
end

"""
mu_P2_per_point_binned_periodic(X, Y, Z, rbin; L=(Lx,Ly,Lz), progress=false, prog_batch=512, updates=100)

Periodic box only. For each point i, bin neighbors by r using `rbin` (Nr+1 edges)
and compute ⟨μ⟩ and ⟨P₂(μ)⟩ per bin. Returns N×(3+2*Nr) Float32:
row i = [x_i, y_i, z_i, ⟨μ⟩₁, ⟨P₂⟩₁, ⟨μ⟩₂, ⟨P₂⟩₂, …].
"""
function mu_P2_per_point_binned_periodic(X, Y, Z, rbin;
                                         L::Tuple=(2000f0,2000f0,2000f0),
                                         progress::Bool=false,
                                         prog_batch::Int=512,
                                         updates::Int=100)

    N  = length(X)
    Nr = length(rbin) - 1

    Xf = Float32.(X); Yf = Float32.(Y); Zf = Float32.(Z)
    rbe = Vector{Float32}(rbin)
    Lx, Ly, Lz = Float32.(L)
    rmax = rbe[end]

    # KDTree wants points as 3×N (columns are points)
    pts = Matrix{Float32}(undef, 3, N)
    @inbounds for i in 1:N
        pts[1,i] = Xf[i]; pts[2,i] = Yf[i]; pts[3,i] = Zf[i]
    end
    tree = KDTree(pts; leafsize=32)

    # Output: (x,y,z) + 2*Nr stats
    out = Array{Float32}(undef, N, 3 + 2*Nr)
    @inbounds for i in 1:N
        out[i,1] = Xf[i]; out[i,2] = Yf[i]; out[i,3] = Zf[i]
    end

    done = progress ? Threads.Atomic{Int}(0) : nothing
    report_every = progress ? max(1, N ÷ max(1, updates)) : 0

    @threads for i in 1:N
        xi = Xf[i]; yi = Yf[i]; zi = Zf[i]

        # IMPORTANT: no kwargs, no StaticArrays
        idxs = inrange(tree, (xi, yi, zi), rmax)

        sμ  = zeros(Float32, Nr)
        sP2 = zeros(Float32, Nr)
        cnt = zeros(Int32,   Nr)

        @inbounds for j in idxs
            j == i && continue
            dx = minimg(Xf[j]-xi, Lx)
            dy = minimg(Yf[j]-yi, Ly)
            dz = minimg(Zf[j]-zi, Lz)

            r2 = muladd(dx, dx, muladd(dy, dy, dz*dz))
            r2 == 0f0 && continue
            r  = sqrt(r2)

            ir = binidx(rbe, r)
            ir == 0 && continue

            μ = abs(dz) / r
            sμ[ir]  += μ
            sP2[ir] += P2(μ)
            cnt[ir] += 1
        end

        base = 3
        @inbounds for k in 1:Nr
            if cnt[k] == 0
                out[i, base + 2k - 1] = NaN32
                out[i, base + 2k]     = NaN32
            else
                invc = 1.0f0 / cnt[k]
                out[i, base + 2k - 1] = sμ[k]  * invc
                out[i, base + 2k]     = sP2[k] * invc
            end
        end

        if progress && (i % prog_batch == 0)
            c = atomic_add!(done, prog_batch)
            (threadid() == 1 && (c % report_every == 0)) && println("Progress: $(min(100, Int(round(100c/N))))%")
        end
    end

    progress && println("Progress: 100%")
    return out
end
