# todo this should be replaced by general broadcasting for <: AbstractGrid
for Grid in (:FullGaussianGrid, :FullClenshawGrid, :FullHEALPixGrid, :FullOctaHEALPixGrid,
            :OctahedralGaussianGrid, :OctahedralClenshawGrid, :HEALPixGrid,
            :OctaHEALPixGrid)
    @eval begin

        function Base.similar(G::$Grid)
            return $Grid{eltype(G)}(undef, G.nlat_half)
        end

        function Base.similar(G::$Grid, ::Type{T}) where T
            return $Grid{T}(undef, G.nlat_half)
        end

        function Base.similar(G::$Grid, nlat_half::Integer)
            return $Grid{eltype(G)}(undef, nlat_half)
        end

        function Base.similar(G::$Grid, ::Type{T}, nlat_half::Integer) where T
            return $Grid{T}(undef, nlat_half)
        end

    end
end