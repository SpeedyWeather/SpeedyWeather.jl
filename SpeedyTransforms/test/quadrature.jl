using Logging
using LinearAlgebra

@testset "Quadrature weights" begin
    HEALPix_grids = (HEALPixGrid, OctaHEALPixGrid, FullHEALPixGrid, FullOctaHEALPixGrid)
    exact_grids = (FullGaussianGrid, OctahedralGaussianGrid, FullClenshawGrid, OctahedralClenshawGrid)

    # the geometric solid angle ΔΩ of ring j, recovered from the transform's fused
    # ΔΩ conj(lon_offset) (the rotation has unit modulus)
    solid_angles(S) = abs.(Array(S.solid_angles_rotated))

    # total solid angle g_j of ring pair j (north + south), at order m; Σ_j g_j = 4π
    function ring_totals(S, m = 1)
        (; nlat, nlons) = S
        ΔΩ = solid_angles(S)
        return [((nlat - j + 1 == j) ? 1 : 2) * nlons[j] * ΔΩ[m, j] for j in 1:S.grid.nlat_half]
    end

    function default_transform(Grid, truncation, NF)
        dealiasing = SpeedyTransforms.default_dealiasing(Grid)
        grid = Grid(SpeedyTransforms.get_nlat_half(truncation, dealiasing))
        return SpectralTransform(Spectrum(truncation), grid; NF)
    end

    @testset "weights are positive and close to equal area" begin
        @testset for Grid in HEALPix_grids
            @testset for truncation in (32, 64)
                S = default_transform(Grid, truncation, Float64)
                ΔΩ = solid_angles(S)
                geometric = RingGrids.get_solid_angles(S.grid)
                @test all(>(0), ΔΩ)
                for j in axes(ΔΩ, 2), m in axes(ΔΩ, 1)
                    @test SpeedyTransforms.QUADRATURE_MIN_RATIO <=
                        ΔΩ[m, j] / geometric[j] <= SpeedyTransforms.QUADRATURE_MAX_RATIO
                end
            end
        end
    end

    @testset "global mean is exactly conserved" begin
        # the l=m=0 exactness condition is Σ_j g_j = 4π over all rings, with g_j the total solid
        # angle of ring pair j; anything else would bias the area-weighted mean of every field
        @testset for Grid in HEALPix_grids
            S = default_transform(Grid, 32, Float64)
            ΔΩ = solid_angles(S)
            (; nlat, nlons) = S
            nlat_half = S.grid.nlat_half
            total = sum(((nlat - j + 1 == j) ? 1 : 2) * nlons[j] * ΔΩ[1, j] for j in 1:nlat_half)
            @test total ≈ 4π
        end
    end

    @testset "no ring is kept at its own Nyquist order" begin
        # a ring at m = nlon/2 carries only one of the two Fourier components, which makes an
        # exact transform unreachable; the HEALPix grids must truncate below it
        @testset for Grid in HEALPix_grids
            S = default_transform(Grid, 32, Float64)
            for j in 1:S.grid.nlat_half
                @test S.mmax_truncation[j] < S.nlons[j] ÷ 2
            end
        end
    end

    @testset "grids with an exact quadrature rule are untouched" begin
        # Gaussian and Clenshaw-Curtis latitudes already carry a quadrature rule, so their weights
        # must stay the grid's geometric solid angles, identical at every order m
        @testset for Grid in exact_grids
            @testset for NF in (Float32, Float64)
                S = default_transform(Grid, 32, NF)
                @test SpeedyTransforms.default_quadrature(Grid) == SpeedyTransforms.EqualAreaQuadrature
                ΔΩ = solid_angles(S)
                geometric = RingGrids.get_solid_angles(S.grid)
                for j in axes(ΔΩ, 2)
                    @test all(ΔΩ[m, j] ≈ NF(geometric[j]) for m in axes(ΔΩ, 1))
                end
            end
        end
    end

    @testset "quadrature schemes" begin
        @testset for Grid in HEALPix_grids
            @testset for truncation in (32, 64)
                nlat_half = SpeedyTransforms.get_nlat_half(truncation, SpeedyTransforms.default_dealiasing(Grid))
                spectrum = Spectrum(truncation)
                geometric = RingGrids.get_solid_angles(Grid(nlat_half))

                # EqualArea reproduces the grid's geometric solid angles exactly, at every order:
                # this is the pre-existing behaviour and must stay bit-identical to it
                S = SpectralTransform(spectrum, Grid(nlat_half); NF = Float64, Quadrature = SpeedyTransforms.EqualAreaQuadrature)
                ΔΩ = solid_angles(S)
                @test all(ΔΩ[m, j] ≈ geometric[j] for j in axes(ΔΩ, 2), m in axes(ΔΩ, 1))

                # Ring weights are one weight per ring, identical at every order …
                S = SpectralTransform(spectrum, Grid(nlat_half); NF = Float64, Quadrature = SpeedyTransforms.RingQuadrature)
                ΔΩ = solid_angles(S)
                @test all(ΔΩ[m, j] ≈ ΔΩ[1, j] for j in axes(ΔΩ, 2), m in axes(ΔΩ, 1))
                @test all(>(0), ΔΩ)
                @test ΔΩ != solid_angles(
                    SpectralTransform(
                        spectrum, Grid(nlat_half); NF = Float64,
                        Quadrature = SpeedyTransforms.EqualAreaQuadrature
                    )
                )

                # … fixed by their defining condition: the quadrature integrates every band-limited
                # function exactly, Σ_j g_j λ_l0(μ_j) = 2√π δ_l0. Only even degrees are meaningful
                # for a north+south ring total; odd ones cancel between the hemispheres.
                λ₀ = Array(S.legendre_polynomials.data)[LowerTriangularArrays.get_lm_range(1, spectrum.lmax - 1), :]
                integrals = (λ₀ * ring_totals(S))[1:2:end]
                @test integrals[1] ≈ 2sqrt(π)                       # ∫Y₀₀ dΩ, i.e. Σ_j g_j = 4π
                @test maximum(abs, integrals[2:end]) < 1.0e-12       # every higher degree integrates to 0

                # the equal-area weights do NOT satisfy that, which is what the scheme fixes
                S₀ = SpectralTransform(spectrum, Grid(nlat_half); NF = Float64, Quadrature = SpeedyTransforms.EqualAreaQuadrature)
                @test maximum(abs, (λ₀ * ring_totals(S₀))[3:2:end]) > 1.0e-5
            end
        end

        # Dropping the ring's Nyquist bin follows the grid, not the scheme, so all three schemes
        # see the same rings and comparing them isolates the weights.
        @testset for Grid in HEALPix_grids
            @test SpeedyTransforms.drop_nyquist(Grid(SpeedyTransforms.get_nlat_half(32, 3.5)))
        end
        @testset for Grid in exact_grids
            @test !SpeedyTransforms.drop_nyquist(Grid(SpeedyTransforms.get_nlat_half(32, 3.0)))
        end
        @testset for Grid in HEALPix_grids
            truncations = [
                SpectralTransform(
                    Spectrum(32), Grid(SpeedyTransforms.get_nlat_half(32, 3.5));
                    NF = Float64, Quadrature = Q
                ).mmax_truncation
                    for Q in (
                        SpeedyTransforms.EqualAreaQuadrature, SpeedyTransforms.RingQuadrature,
                        SpeedyTransforms.PerOrderQuadrature, SpeedyTransforms.ContractiveQuadrature,
                    )
            ]
            @test allequal(truncations)
        end
    end

    @testset "ContractiveQuadrature makes analysis non-expansive" begin
        # In the physical (solid-angle) norm on grid space the analysis operator at order m is
        # A = Λ diag(g) diag(g⁰)^(-1/2), block diagonal in the parity of l-m. σmax(A) ≤ 1 means no
        # grid field — in band or aliased above the truncation — can be analysed into more spectral
        # energy than it had, which is the property PerOrderQuadrature gives up to gain exactness.
        function analysis_norm(S)
            (; nlat, nlons, mmax_truncation) = S
            (; lmax, mmax) = S.spectrum
            nlat_half = S.grid.nlat_half
            λ = Array(S.legendre_polynomials.data)
            ΔΩ = solid_angles(S)
            multiplicity(j) = (nlat - j + 1 == j) ? 1 : 2
            g⁰ = [multiplicity(j) * nlons[j] * Float64(RingGrids.get_solid_angles(S.grid)[j]) for j in 1:nlat_half]
            σmax = 0.0
            for m in 1:mmax
                rings = [j for j in 1:nlat_half if mmax_truncation[j] + 1 >= m]
                isempty(rings) && continue
                Λ = Float64.(λ[LowerTriangularArrays.get_lm_range(m, lmax - 1), rings])
                g = [multiplicity(j) * nlons[j] * Float64(ΔΩ[m, j]) for j in rings]
                A = Λ * LinearAlgebra.Diagonal(g ./ sqrt.(g⁰[rings]))
                for parity in (0, 1)
                    block = @view A[(1 + parity):2:end, :]
                    size(block, 1) == 0 && continue
                    σmax = max(σmax, LinearAlgebra.opnorm(block))
                end
            end
            return σmax
        end

        @testset for Grid in HEALPix_grids
            @testset for truncation in (32, 64)
                nlat_half = SpeedyTransforms.get_nlat_half(truncation, SpeedyTransforms.default_dealiasing(Grid))
                spectrum = Spectrum(truncation)
                transform(Q) = SpectralTransform(spectrum, Grid(nlat_half); NF = Float64, Quadrature = Q)

                # the property the scheme exists for, to roundoff
                @test analysis_norm(transform(SpeedyTransforms.ContractiveQuadrature)) <= 1 + 1.0e-12

                # it is a strict improvement on both other HEALPix schemes, which amplify
                @test analysis_norm(transform(SpeedyTransforms.EqualAreaQuadrature)) > 1 + 1.0e-6
                @test analysis_norm(transform(SpeedyTransforms.PerOrderQuadrature)) > 1 + 1.0e-3

                # weights stay positive and close to equal area. They no longer only shrink: the
                # fit moves them either way, by up to ~9% at the sanctioned dealiasing, and the
                # shrink that enforces ‖A‖ ≤ 1 is applied on top of that.
                geometric = RingGrids.get_solid_angles(Grid(nlat_half))
                ΔΩ = solid_angles(transform(SpeedyTransforms.ContractiveQuadrature))
                @test all(>(0), ΔΩ)
                for j in axes(ΔΩ, 2), m in axes(ΔΩ, 1)
                    @test 0.9 <= ΔΩ[m, j] / geometric[j] <= 1.2
                end
            end
        end

        # The point of fitting rather than uniformly rescaling: at the sanctioned dealiasing the
        # round trip must beat EqualAreaQuadrature, which the plain scalar shrink does not. Compared
        # in grid space, as in the "Transform roundtrip accuracy" testset, because raw random
        # coefficients carry components no real grid field can represent.
        @testset "more accurate than equal area" begin
            @testset for Grid in (HEALPixGrid, OctaHEALPixGrid)
                @testset for truncation in (32, 64)
                    nlat_half = SpeedyTransforms.get_nlat_half(
                        truncation, SpeedyTransforms.default_dealiasing(Grid)
                    )
                    spectrum = Spectrum(truncation)
                    spec = randn(Complex{Float64}, spectrum)

                    function roundtrip_error(Quadrature)
                        S = SpectralTransform(spectrum, Grid(nlat_half); NF = Float64, Quadrature)
                        field = transform(spec, S)
                        roundtrip = transform(transform(field, S), S)
                        data, reference = Array(roundtrip.data), Array(field.data)
                        return norm(data - reference) / norm(reference)
                    end

                    equal_area = roundtrip_error(SpeedyTransforms.EqualAreaQuadrature)
                    contractive = roundtrip_error(SpeedyTransforms.ContractiveQuadrature)
                    @test contractive < equal_area
                end
            end
        end

        # grids that already carry an exact quadrature rule are already non-expansive, so the scheme
        # must leave their geometric solid angles untouched
        @testset for Grid in exact_grids
            nlat_half = SpeedyTransforms.get_nlat_half(32, SpeedyTransforms.default_dealiasing(Grid))
            S = SpectralTransform(
                Spectrum(32), Grid(nlat_half);
                NF = Float64, Quadrature = SpeedyTransforms.ContractiveQuadrature
            )
            geometric = RingGrids.get_solid_angles(Grid(nlat_half))
            ΔΩ = solid_angles(S)
            @test all(ΔΩ[m, j] ≈ geometric[j] for j in axes(ΔΩ, 2), m in axes(ΔΩ, 1))
        end
    end

    @testset "an under-resolved grid warns and falls back" begin
        # dealiasing 2 leaves whole orders m uncovered by any ring on HEALPix; no weighting can
        # fix that, so the fit must decline rather than return unusable weights
        grid = HEALPixGrid(SpeedyTransforms.get_nlat_half(64, 2))
        S = @test_logs (:warn,) match_mode = :any SpectralTransform(Spectrum(64), grid; NF = Float64)
        ΔΩ = solid_angles(S)
        geometric = RingGrids.get_solid_angles(S.grid)
        @test all(>(0), ΔΩ)
        @test all(
            ΔΩ[m, j] <= SpeedyTransforms.QUADRATURE_MAX_RATIO * geometric[j]
                for j in axes(ΔΩ, 2), m in axes(ΔΩ, 1)
        )
    end

    @testset "default dealiasing" begin
        @testset for Grid in HEALPix_grids
            @test SpeedyTransforms.default_dealiasing(Grid) == 3.5
            # dealiasing 3.5 puts the target resolutions on grids HEALPixGrid can represent, with
            # a spectral slack of nlat_half/9 — enough for a well-conditioned exact fit on both
            # HEALPix families, see docs/dev/2026-08/healpix-quadrature-exactness.md
            for truncation in (32, 64, 128, 256)
                nlat_half = SpeedyTransforms.get_nlat_half(truncation, 3.5)
                @test iseven(nlat_half)                     # HEALPixGrid requires even nlat_half
                @test nlat_half - truncation >= cld(nlat_half, 9)
            end
        end
        @testset for Grid in exact_grids
            @test SpeedyTransforms.default_dealiasing(Grid) < 3.5
        end
    end

    @testset "an inexact fit is reported even when the weights look reasonable" begin
        # nlat_half ≥ truncation here, so the system is solvable and the ring-count fallback does
        # not trigger, but it is too ill-conditioned to reach exactness. Only the residual check
        # catches that — the weight bounds alone would let it pass silently.
        logger = Test.TestLogger()
        S = with_logger(logger) do
            SpectralTransform(Spectrum(71), OctaHEALPixGrid(72); NF = Float64)
        end
        @test any(r -> r.level >= Logging.Warn, logger.logs)
        @test all(>(0), solid_angles(S))        # still usable weights, just not exact
    end

    @testset "every order is fitted at the default dealiasing" begin
        # the fit warns when it has to fall back on any order, so a silent construction is the
        # direct check that the default dealiasing leaves enough slack at higher truncations too
        @testset for Grid in (HEALPixGrid, OctaHEALPixGrid)
            @test_logs min_level = Logging.Warn default_transform(Grid, 128, Float64)
        end
    end
end
