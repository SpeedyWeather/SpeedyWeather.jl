import Pkg

@testset "Asset bindings by source and version" begin
    asset_path = joinpath("data", "boundary_conditions", "orography.nc")
    source = "https://example.com/assets"
    binding_v1 = RingGrids.asset_binding_name(asset_path; assets_url = source, version = v"1")
    binding_v2 = RingGrids.asset_binding_name(asset_path; assets_url = source, version = v"2")
    binding_branch = RingGrids.asset_binding_name(asset_path; assets_url = source, version = "main")
    binding_other_repo = RingGrids.asset_binding_name(
        asset_path; assets_url = "https://example.com/other-assets", version = v"1"
    )
    binding_other_path = RingGrids.asset_binding_name(
        joinpath("data", "other", "orography.nc"); assets_url = source, version = v"1"
    )
    binding_tag = RingGrids.asset_binding_name(asset_path; assets_url = source, version = v"1.0.0")
    binding_named_branch = RingGrids.asset_binding_name(asset_path; assets_url = source, version = "v1.0.0")

    @test binding_v1 == joinpath(source, "tags/v1.0.0", asset_path)
    @test binding_branch == joinpath(source, "heads/main", asset_path)
    @test binding_v1 == binding_tag
    @test binding_tag != binding_named_branch
    @test length(unique((binding_v1, binding_v2, binding_branch, binding_other_repo, binding_other_path))) == 5
    @test RingGrids.asset_binding_name(asset_path) == RingGrids.asset_binding_name(
        asset_path; assets_url = RingGrids.ASSETS_URL, version = RingGrids.DEFAULT_ASSETS_VERSION
    )

    mktempdir() do dir
        artifact_toml = joinpath(dir, "Artifacts.toml")
        hash_v1 = Base.SHA1("0123456789abcdef0123456789abcdef01234567")
        hash_v2 = Base.SHA1("abcdef0123456789abcdef0123456789abcdef01")

        # A legacy filename binding cannot be reused for an unknown version.
        Pkg.Artifacts.bind_artifact!(artifact_toml, basename(asset_path), hash_v1)
        @test Pkg.Artifacts.artifact_hash(binding_v1, artifact_toml) === nothing

        # Binding another version leaves earlier versions available.
        Pkg.Artifacts.bind_artifact!(artifact_toml, binding_v1, hash_v1)
        @test Pkg.Artifacts.artifact_hash(binding_v2, artifact_toml) === nothing
        Pkg.Artifacts.bind_artifact!(artifact_toml, binding_v2, hash_v2)
        @test Pkg.Artifacts.artifact_hash(binding_v1, artifact_toml) == hash_v1
        @test Pkg.Artifacts.artifact_hash(binding_v2, artifact_toml) == hash_v2

        for binding in (binding_branch, binding_other_repo, binding_other_path)
            @test Pkg.Artifacts.artifact_hash(binding, artifact_toml) === nothing
            Pkg.Artifacts.bind_artifact!(artifact_toml, binding, hash_v2)
            @test Pkg.Artifacts.artifact_hash(binding, artifact_toml) == hash_v2
            @test Pkg.Artifacts.artifact_hash(binding_v1, artifact_toml) == hash_v1
        end

        # Looking up a previous version reads its saved binding without changes.
        original = read(artifact_toml, String)
        @test Pkg.Artifacts.artifact_hash(binding_tag, artifact_toml) == hash_v1
        @test read(artifact_toml, String) == original
        @test Pkg.Artifacts.artifact_hash(basename(asset_path), artifact_toml) == hash_v1

        # Refreshing one branch leaves other versions and repositories cached.
        Pkg.Artifacts.unbind_artifact!(artifact_toml, binding_branch)
        @test Pkg.Artifacts.artifact_hash(binding_branch, artifact_toml) === nothing
        @test Pkg.Artifacts.artifact_hash(binding_v1, artifact_toml) == hash_v1
        @test Pkg.Artifacts.artifact_hash(binding_v2, artifact_toml) == hash_v2
        @test Pkg.Artifacts.artifact_hash(binding_other_repo, artifact_toml) == hash_v2
    end
end
