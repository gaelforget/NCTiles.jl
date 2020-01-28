
"""
    findidx(A,val)

Helper function for getting the indices for tiles. A is a gcmfaces struct, val is a
    numeric value to get the indices of. Returns a Dict of the indices for each face.
    Not currently exported. Maybe move to MeshArrays?
"""
function findidx(A,val)
    tileidx = Dict()
    if isa(A,MeshArrays.gcmfaces)
        nfaces = A.nFaces
    elseif isa(A,MeshArray)
        nfaces = A.grid.nFaces
    end
    for i = 1:nfaces
        idx = findall(x->x==val,A.f[i])
        tileidx[i] = [x for x in idx]
    end
    return tileidx
end

"""
    gettile(fldvals,tileinfo,tilesize,tilenum::Int)

Helper function for retrieving a tile from a gcmfaces struct as a numeric Array. Not
    currently exported.
"""
function gettile(fldvals,tileinfo,tilesize,tilenum::Int)
    tilidx = findidx(tileinfo["tileNo"],tilenum)

    if isa(fldvals,MeshArrays.gcmfaces)
        is3D = length(size(fldvals)) == 3
        if is3D; n3 = size(fldvals)[3]; end
    else
        is3D = length(size(fldvals)) == 2 && size(fldvals)[2] > 1
        if is3D; n3 = size(fldvals)[2]; end
    end

    if is3D # has depth
        tilfld = Array{Any,3}(nothing,0,tilesize[2],n3)
    else
        tilfld = Array{Any,2}(nothing,0,tilesize[2])
    end
    if isa(fldvals,MeshArrays.gcmfaces)
        nfaces = fldvals.nFaces
    elseif isa(fldvals,MeshArray)
        nfaces = fldvals.grid.nFaces
    end
    for iF = 1:nfaces
        if ~isempty(tilidx[iF])
            imin = minimum(tilidx[iF])[1]; imax = maximum(tilidx[iF])[1]
            jmin = minimum(tilidx[iF])[2]; jmax = maximum(tilidx[iF])[2]
            if is3D
                if isa(fldvals,MeshArray)
                    tilfld = [tilfld; cat([fldvals.f[iF,d][imin:imax,jmin:jmax] for d in 1:n3]...,dims=3)]
                else
                    tilfld = [tilfld; fldvals.f[iF][imin:imax,jmin:jmax,:]]
                end
            else
                tilfld = [tilfld; fldvals.f[iF][imin:imax,jmin:jmax]]
            end
        end
    end

    if is3D
        tilfld = reshape(tilfld,tilesize[1],tilesize[2],:)
    else
        tilfld = reshape(tilfld,tilesize[1],tilesize[2])
    end

    return tilfld
end

"""
    gettiles(fldvals,tilenum::Int)

Helper function for retrieving a tile from a gcmfaces struct as a numeric Array along
    with associated latitude and longitude. Not currently exported.
"""
function gettiles(tilfld,tilenum::Int)
    tilesize = tilfld.tilesize

    tilfld = gettile(tilfld.vals,tilfld.tileinfo,tilesize,tilenum)
    tillat = gettile(tilfld.tileinfo["XC"],tilfld.tileinfo,tilfld.tilesize,tilenum)
    tillon = gettile(tilfld.tileinfo["YC"],tilfld.tileinfo,tilfld.tilesize,tilenum)

    return tilfld,tillat,tillon
end
