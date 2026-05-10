using CUDA

# --- CPU ---
A_cpu = zeros(Float32, 5, 6)   # 5 grid points, 6 layers

# view() on CPU returns a SubArray sharing the same backing memory
v_cpu = view(A_cpu, :, 1:3)
println("CPU view type:     ", typeof(v_cpu))       # SubArray{...}
println("CPU Base.parent:   ", Base.parent(v_cpu) === A_cpu)  # true

# Write through parent → visible in view (and vice versa)
A_cpu[:, 1:3] .= 7f0
println("CPU parent→view:   ", all(v_cpu .== 7f0))  # true
v_cpu .= 0f0
println("CPU view→parent:   ", all(A_cpu[:, 1:3] .== 0f0))  # true

# --- GPU ---
A_gpu = CUDA.zeros(Float32, 5, 6)

# view() on a CuArray returns a NEW CuArray over the same GPU buffer (offset pointer)
# It is NOT a SubArray, so Base.parent(v_gpu) !== A_gpu
v_gpu = view(A_gpu, :, 1:3)
println("\nGPU view type:     ", typeof(v_gpu))         # CuArray{...} — NOT SubArray
println("GPU Base.parent:   ", Base.parent(v_gpu) === A_gpu)  # false ← this is what trips the test

# But the memory IS shared — writes go both ways, just like SubArray
A_gpu[:, 1:3] .= 7f0
println("GPU parent→view:   ", all(Array(v_gpu) .== 7f0))   # true  ← memory is shared
v_gpu .= 0f0
println("GPU view→parent:   ", all(Array(A_gpu[:, 1:3]) .== 0f0))  # true

# The correct portable check for "is this memory shared?" is behavioral, not structural:
println("\nPortable memory-sharing check: fill parent, see value through view")
A_gpu .= 0f0
A_gpu[:, 1:3] .= 42f0
println("GPU shared memory: ", all(Array(v_gpu) .== 42f0))   # true — same buffer
