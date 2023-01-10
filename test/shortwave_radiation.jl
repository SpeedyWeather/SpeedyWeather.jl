@testset "Parametrization: shortwave radiation" begin
    @testset "solar!" begin
        @testset for NF in (Float32, Float64)
            topsr = [0., 283.06, 352.56537, 378.8209880964894]
            latd = [-89., 0., 45., 89.]
            @testset for i in 1:4
                _, diagn, model = SpeedyWeather.initialize_speedy(NF, PrimitiveEquation, nlev = 4)
                nlev = diagn.nlev
                nband = model.parameters.nband

                column = ColumnVariables{NF}(nlev=nlev, nband=nband)
                column.tyear = 0.5
                column.csol = 1000.
                column.latd = latd[i]

                SpeedyWeather.solar!(column)

                @test column.topsr ≈ topsr[i]
            end
        end
    end
    @testset "sol_oz!" begin
        @testset for NF in (Float32, Float64)
            _, diagn, model = SpeedyWeather.initialize_speedy(NF, PrimitiveEquation, nlev = 4)
            nlev = diagn.nlev
            nband = model.parameters.nband

            column = ColumnVariables{NF}(nlev=nlev, nband=nband)
            column.tyear = 0.5
            column.latd = 89.

            SpeedyWeather.sol_oz!(column, model)

            @test column.topsr ≈ 518.2271117159975
            @test column.fsol ≈ 518.2271117159975
            @test column.ozupp ≈ 6.99611021887365
            @test column.ozone ≈ 15.666684
            @test column.zenit ≈ 1.3500085311452612
            @test column.stratz ≈ 0.
        end
    end
    @testset "cloud!" begin
        @testset for NF in (Float32, Float64)
            _, diagn, model = SpeedyWeather.initialize_speedy(NF, PrimitiveEquation, nlev = 4)
            nlev = diagn.nlev
            nband = model.parameters.nband

            column = ColumnVariables{NF}(nlev=nlev, nband=nband)
            column.humid .= [0., 0.03124937, 0.9748285, 6.7846994] # g/kg
            column.rel_hum = [0.0, 1., 1., 0.8964701087948754]
            column.grad_dry_static_energy = 0.4255393541723314
            column.precip_convection = 1
            column.precip_large_scale = 1
            column.cloud_top = 1
            column.fmask = 1.

            SpeedyWeather.cloud!(column, model)

            @test column.icltop ≈ 1.
            @test column.cloudc ≈ 1.
            @test column.clstr ≈ 0.13447051631923132
            @test column.qcloud ≈ 0.9748285
        end
    end
    @testset "radsw!" begin
        @testset for NF in (Float32, Float64)
            _, diagn, model = SpeedyWeather.initialize_speedy(NF, PrimitiveEquation, nlev = 4)
            nlev = diagn.nlev
            nband = model.parameters.nband

            column = ColumnVariables{NF}(nlev=nlev, nband=nband)
            column.norm_pres = 1.
            column.humid .= [0., 0.03124937, 0.9748285, 6.7846994] # g/kg
            column.albsfc = 0.5
            # Same as returned from sol_oz and cloud
            column.icltop = 1
            column.cloudc = 0.90186596
            column.clstr = 0.1485
            column.ozupp = 6.99611021887365
            column.ozone = 15.666684
            column.zenit = 1.3500085311452612
            column.stratz = 0.
            column.fsol = 518.2271117159975
            column.qcloud = 0.033334

            SpeedyWeather.radsw!(column, model)

            @test column.ssrd ≈ 385.7997293028816
            @test column.ssr ≈ 192.8998646514408
            @test column.tsr ≈ 315.10016
            @test column.tau2 ≈ [0.9526258427031802 0.3788324531779398 1.0 1.0;
                                 0.9286716429916902 0.2276374390970355 0.994618802250959 0.6801722651144487;
                                 0.02148265489609933 0.12596241035550795 0.7900788064253661 4.906182803130207e-8;
                                 0.9287847242363858 0.22819245387064319 0.31050206404533354 5.234705430612321e-37]
            @test column.tend_t_rsw ≈ [11.947794979452055, 27.611640775148402, 42.565299927060124, 40.46336031165737]
        end
    end
    @testset "shortwave_radiation!" begin
        @testset for NF in (Float32, Float64)
            _, diagn, model = SpeedyWeather.initialize_speedy(NF, PrimitiveEquation, nlev = 4)
            nlev = diagn.nlev
            nband = model.parameters.nband

            column = ColumnVariables{NF}(nlev=nlev, nband=nband)

            # 1. Set variables for sol_oz
            column.tyear = 0.5
            column.latd = 89.

            # 2. Compute sat_vap_pres and dry_static_energy
            # and set remainig varables for cloud
            column.pres = 1. # nomalised
            column.temp .= [208.40541, 219.8126 , 249.25502, 276.14264]
            column.humid .= [0., 0.03124937, 0.9748285, 6.7846994] # g/kg
            column.geopot .= [241932.19, 117422.14, 54618.79, 7626.5884]
            SpeedyWeather.get_thermodynamics!(column, model)

            # column.sat_vap_pres = [0.005265206274688095, 0.02494611200040165, 0.7012750147882735, 7.56823828640608]
            # column.dry_static_energy = [451171.22164, 338113.9904, 304870.83008, 284873.79896]
            # column.rel_hum = [0.0, 1.2526749659224194, 1.3900801817305823, 0.8964701087948754]
            # column.grad_dry_static_energy = 0.4255393541723314

            column.precip_convection = 1
            column.precip_large_scale = 1
            column.cloud_top = 1
            column.fmask = 1.

            # Set variables for radsw
            # FIXME: pres is overloaded and will need to be
            # fixed in other functions
            column.norm_pres = column.pres
            column.albsfc = 0.5

            SpeedyWeather.shortwave_radiation!(column, model)

            # Results need to be close to those computed by radsw but
            # are not the same as precision is maintained between functions
            @test isapprox(column.ssrd, 385.7997293028816, rtol=0.1)
            @test isapprox(column.ssr, 192.8998646514408, rtol=0.1)
            @test isapprox(column.tsr, 315.10016, rtol=0.1)
            @test isapprox(column.tau2, [
                0.9526258427031802 0.3788324531779398 1.0 1.0;
                0.9286716429916902 0.2276374390970355 0.994618802250959 0.6801722651144487;
                0.02148265489609933 0.12596241035550795 0.7900788064253661 4.906182803130207e-8;
                0.9287847242363858 0.22819245387064319 0.31050206404533354 5.234705430612321e-37
                ], rtol=0.1)
            @test isapprox(column.tend_t_rsw, [
                11.947794979452055, 27.611640775148402, 42.565299927060124, 40.46336031165737
                ], rtol=0.1)
            # Range is generally order 10^-5
            @test column.temp_tend ≈ [7.189150408440864e-6, 1.2141337616807728e-5, 1.3196339560843619e-5, 1.5994200814647754e-5]
        end
    end
end
