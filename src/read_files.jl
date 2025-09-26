using FITSIO

"""
    read_xyz_from_dir(dir::AbstractString) -> Matrix{Float32}
    read_xyz(files::Vector{String}) -> Matrix{Float32}

Read EZmock `.fits.gz` files, extract X,Y,Z columns (case-insensitive),
concatenate, and return as N×3 Float32 matrix.
"""
function read_xyz_from_dir(dir::AbstractString)
    files = filter(f -> occursin(r"^EZmock.*\.fits\.gz$"i, basename(f)),
                   readdir(dir; join=true))
    return read_xyz(files)
end

function read_xyz(files::Vector{String})
    Xall = Float32[]; Yall = Float32[]; Zall = Float32[]

    for file in files
        FITS(file) do f
            # find a table HDU that has columns X,Y,Z (case-insensitive)
            found = false
            for hdu in f
                if hdu isa FITSIO.TableHDU
                    names = FITSIO.colnames(hdu)  # Vector{String}
                    up = uppercase.(names)
                    ix = findfirst(==("X"), up)
                    iy = findfirst(==("Y"), up)
                    iz = findfirst(==("Z"), up)
                    if !(ix === nothing || iy === nothing || iz === nothing)
                        # use the original-cased names when reading
                        X = read(hdu, names[ix])
                        Y = read(hdu, names[iy])
                        Z = read(hdu, names[iz])
                        append!(Xall, Float32.(X))
                        append!(Yall, Float32.(Y))
                        append!(Zall, Float32.(Z))
                        found = true
                        break
                    end
                end
            end
            found || error("No table HDU with X,Y,Z in file: $file")
        end
    end

    return hcat(Xall, Yall, Zall)
end
