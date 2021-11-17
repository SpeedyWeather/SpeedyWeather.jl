"""
Perform a leapfrog timestep using a RAW filter given field at t=0 and derivatives
"""
function leapfrog( f::Array{Complex{T},3},           # input/output field of dimension m * n * 2. The "2" is required for timesetepping 
                   df::Array{Complex{T},2},          # derivatives of field f of dimension m * n
                   LF:: LeapfrogSettings{T}          # struct holding integration parameters
                   ) where {T<:AbstractFloat}

    #Get the dimensions of f             
    mx,nx,2 = size(f)

    #Check the dimensions of the derivative field match OK
    @boundscheck (mx,nx) == size(df)   || throw(BoundsError())

    #Get the timestep, Williams parameter and RA parameter
    @unpack dt, wil, eps = LF.settings


    #ix = 96, iy = 24
    #if (ix == iy*4) then
       # call trunct(fdt) # vor = vor * trfilt. Some kind of rhomboidal truncation
    #end if

    #The actual leap frog with the Robert filter
    fnew = f(:,:,1) + dt*df
    f(:,:,1) = f(:,:,2) + wil*eps*(f(:,:,1) - 2*output(:,:,2) + fnew)

    #Williams' innovation to the filter
    f(:,:,2) = fnew - (1.0 - wil)*eps*(f(:,:,1) - 2.0*f(:,:,2) + fnew)

end
