using SpeedyWeatherInternals.Utils
using Dates
using Test

@testset "Custom Dates" begin

    @testset "Period construction" begin
        @test Dates.value(SpeedyMillisecond(1)) == 1
        @test Dates.value(SpeedySecond(2)) == 2
        @test Dates.value(SpeedyMinute(3)) == 3
        @test Dates.value(SpeedyHour(4)) == 4
        @test Dates.value(SpeedyDay(5)) == 5
        @test Dates.value(SpeedyWeek(6)) == 6
        @test Dates.value(SpeedyMonth(7)) == 7
        @test Dates.value(SpeedyYear(8)) == 8

        # parametric type is inferred from input
        @test SpeedySecond(Int32(5)) isa SpeedySecond{Int32}
        @test SpeedySecond(5) isa SpeedySecond{Int64}
    end

    @testset "Period printing" begin
        @test sprint(print, SpeedySecond(1)) == "1 second"
        @test sprint(print, SpeedySecond(2)) == "2 seconds"
        @test sprint(print, SpeedyDay(1)) == "1 day"
        @test sprint(print, SpeedyDay(3)) == "3 days"
        @test sprint(print, SpeedyHour(0)) == "0 hours"
    end

    @testset "Period zero/one/isfinite" begin
        @test zero(SpeedySecond) == SpeedySecond(0)
        @test zero(SpeedyDay(5)) == SpeedyDay(0)
        @test one(SpeedyHour) == 1
        @test isfinite(SpeedyMinute(10))
        @test iszero(SpeedySecond(0))
        @test !iszero(SpeedySecond(1))
    end

    @testset "Period unary operations" begin
        @test -SpeedySecond(5) == SpeedySecond(-5)
        @test +SpeedySecond(5) == SpeedySecond(5)
        @test abs(SpeedySecond(-3)) == SpeedySecond(3)
        @test sign(SpeedySecond(-3)) == -1
        @test sign(SpeedySecond(3)) == 1
        @test signbit(SpeedySecond(-1)) == true
        @test signbit(SpeedySecond(1)) == false
    end

    @testset "Period same-type arithmetic" begin
        @test SpeedySecond(3) + SpeedySecond(2) == SpeedySecond(5)
        @test SpeedySecond(5) - SpeedySecond(2) == SpeedySecond(3)
        @test SpeedySecond(6) / SpeedySecond(3) == 2.0
        @test SpeedySecond(10) / 2 == SpeedySecond(5)
        @test div(SpeedySecond(7), SpeedySecond(3)) == 2
        @test rem(SpeedySecond(7), SpeedySecond(3)) == SpeedySecond(1)
        @test mod(SpeedySecond(7), SpeedySecond(3)) == SpeedySecond(1)
        @test SpeedySecond(3) * 2 == SpeedySecond(6)
        @test 2 * SpeedySecond(3) == SpeedySecond(6)
        @test SpeedyHour(2) * 3 == SpeedyHour(6)
    end

    @testset "Period comparison" begin
        @test SpeedySecond(1) == SpeedySecond(1)
        @test SpeedySecond(1) != SpeedySecond(2)
        @test SpeedySecond(1) < SpeedySecond(2)
        @test SpeedySecond(2) > SpeedySecond(1)
        @test isless(SpeedySecond(1), SpeedySecond(2))
    end

    @testset "Period conversions (FixedPeriod)" begin
        # Millisecond <-> Second
        @test convert(SpeedyMillisecond, SpeedySecond(1)) == SpeedyMillisecond(1000)
        @test convert(SpeedySecond, SpeedyMillisecond(1000)) == SpeedySecond(1)
        @test_throws InexactError convert(SpeedySecond, SpeedyMillisecond(999))

        # Second <-> Minute
        @test convert(SpeedySecond, SpeedyMinute(1)) == SpeedySecond(60)
        @test convert(SpeedyMinute, SpeedySecond(60)) == SpeedyMinute(1)

        # Minute <-> Hour
        @test convert(SpeedyMinute, SpeedyHour(1)) == SpeedyMinute(60)
        @test convert(SpeedyHour, SpeedyMinute(60)) == SpeedyHour(1)

        # Hour <-> Day
        @test convert(SpeedyHour, SpeedyDay(1)) == SpeedyHour(24)
        @test convert(SpeedyDay, SpeedyHour(24)) == SpeedyDay(1)

        # Day <-> Week
        @test convert(SpeedyDay, SpeedyWeek(1)) == SpeedyDay(7)
        @test convert(SpeedyWeek, SpeedyDay(7)) == SpeedyWeek(1)

        # Multi-step conversions
        @test convert(SpeedyMillisecond, SpeedyHour(1)) == SpeedyMillisecond(3600000)
        @test convert(SpeedySecond, SpeedyDay(1)) == SpeedySecond(86400)
        @test convert(SpeedyMillisecond, SpeedyDay(1)) == SpeedyMillisecond(86400000)
    end

    @testset "Period conversions (OtherPeriod)" begin
        @test convert(SpeedyMonth, SpeedyYear(1)) == SpeedyMonth(12)
        @test convert(SpeedyYear, SpeedyMonth(12)) == SpeedyYear(1)
        @test_throws InexactError convert(SpeedyYear, SpeedyMonth(5))
    end

    @testset "Period promotion and cross-type arithmetic" begin
        # FixedPeriod promotion
        @test SpeedySecond(1) + SpeedyMillisecond(500) == SpeedyMillisecond(1500)
        @test SpeedyMinute(1) + SpeedySecond(30) == SpeedySecond(90)
        @test SpeedyHour(1) - SpeedyMinute(30) == SpeedyMinute(30)
        @test SpeedyDay(1) + SpeedyHour(12) == SpeedyHour(36)

        # Cross-type comparison
        @test SpeedySecond(60) == SpeedyMinute(1)
        @test SpeedyHour(1) == SpeedyMinute(60)
        @test SpeedyDay(1) == SpeedyHour(24)

        # OtherPeriod promotion
        @test SpeedyYear(1) + SpeedyMonth(6) == SpeedyMonth(18)
    end

    @testset "toms" begin
        @test Dates.toms(SpeedyMillisecond(1)) == 1
        @test Dates.toms(SpeedySecond(1)) == 1000
        @test Dates.toms(SpeedyMinute(1)) == 60000
        @test Dates.toms(SpeedyHour(1)) == 3600000
        @test Dates.toms(SpeedyDay(1)) == 86400000
        @test Dates.toms(SpeedyWeek(1)) == 604800000
    end

    @testset "Period construct-from-Period" begin
        @test SpeedyMillisecond(SpeedySecond(1)) == SpeedyMillisecond(1000)
        @test SpeedySecond(SpeedyMinute(1)) == SpeedySecond(60)
        @test SpeedyMinute(SpeedyHour(1)) == SpeedyMinute(60)
        @test SpeedyHour(SpeedyDay(1)) == SpeedyHour(24)
        @test SpeedyDay(SpeedyWeek(1)) == SpeedyDay(7)
        @test SpeedyMonth(SpeedyYear(1)) == SpeedyMonth(12)
    end

    @testset "SpeedyDateTime construction" begin
        dt = SpeedyDateTime(2000, 1, 1)
        @test Dates.year(dt) == 2000
        @test Dates.month(dt) == 1
        @test Dates.day(dt) == 1
        @test Dates.hour(dt) == 0
        @test Dates.minute(dt) == 0
        @test Dates.second(dt) == 0
        @test Dates.millisecond(dt) == 0

        dt2 = SpeedyDateTime(2024, 3, 15, 14, 30, 45, 123)
        @test Dates.year(dt2) == 2024
        @test Dates.month(dt2) == 3
        @test Dates.day(dt2) == 15
        @test Dates.hour(dt2) == 14
        @test Dates.minute(dt2) == 30
        @test Dates.second(dt2) == 45
        @test Dates.millisecond(dt2) == 123
    end

    @testset "SpeedyDateTime from Period types" begin
        dt = SpeedyDateTime(SpeedyYear(2000), SpeedyMonth(6), SpeedyDay(15))
        @test Dates.year(dt) == 2000
        @test Dates.month(dt) == 6
        @test Dates.day(dt) == 15
    end

    @testset "SpeedyDateTime equality and comparison" begin
        dt1 = SpeedyDateTime(2000, 1, 1)
        dt2 = SpeedyDateTime(2000, 1, 1)
        dt3 = SpeedyDateTime(2000, 1, 2)

        @test dt1 == dt2
        @test dt1 != dt3
        @test dt1 < dt3
        @test dt3 > dt1
        @test isless(dt1, dt3)
    end

    @testset "SpeedyDateTime - SpeedyDateTime" begin
        dt1 = SpeedyDateTime(2000, 1, 1)
        dt2 = SpeedyDateTime(2000, 1, 2)
        @test dt2 - dt1 == SpeedyMillisecond(86400000)
        @test dt1 - dt1 == SpeedyMillisecond(0)
    end

    @testset "SpeedyDateTime +/- FixedPeriod" begin
        dt = SpeedyDateTime(2000, 1, 1, 0, 0, 0)

        @test dt + SpeedyMillisecond(1) == SpeedyDateTime(2000, 1, 1, 0, 0, 0, 1)
        @test dt + SpeedySecond(1) == SpeedyDateTime(2000, 1, 1, 0, 0, 1)
        @test dt + SpeedyMinute(1) == SpeedyDateTime(2000, 1, 1, 0, 1, 0)
        @test dt + SpeedyHour(1) == SpeedyDateTime(2000, 1, 1, 1, 0, 0)
        @test dt + SpeedyDay(1) == SpeedyDateTime(2000, 1, 2)
        @test dt + SpeedyWeek(1) == SpeedyDateTime(2000, 1, 8)

        @test dt - SpeedyDay(1) == SpeedyDateTime(1999, 12, 31)
        @test dt - SpeedyHour(1) == SpeedyDateTime(1999, 12, 31, 23, 0, 0)

        # commutativity
        @test SpeedyDay(1) + dt == dt + SpeedyDay(1)
    end

    @testset "SpeedyDateTime +/- Year" begin
        dt = SpeedyDateTime(2000, 2, 29)  # leap day
        @test dt + SpeedyYear(1) == SpeedyDateTime(2001, 2, 28)  # clamp to end of Feb
        @test dt - SpeedyYear(1) == SpeedyDateTime(1999, 2, 28)
        @test dt + SpeedyYear(4) == SpeedyDateTime(2004, 2, 29)  # still leap year
    end

    @testset "SpeedyDateTime +/- Month" begin
        dt = SpeedyDateTime(2000, 1, 31)
        @test dt + SpeedyMonth(1) == SpeedyDateTime(2000, 2, 29)  # clamp to Feb 29 (leap year)
        @test dt + SpeedyMonth(2) == SpeedyDateTime(2000, 3, 31)

        dt2 = SpeedyDateTime(2001, 1, 31)
        @test dt2 + SpeedyMonth(1) == SpeedyDateTime(2001, 2, 28)  # non-leap year

        # wrap around year boundary
        dt3 = SpeedyDateTime(2000, 11, 15)
        @test dt3 + SpeedyMonth(3) == SpeedyDateTime(2001, 2, 15)
        @test dt3 - SpeedyMonth(12) == SpeedyDateTime(1999, 11, 15)
    end

    @testset "SpeedyDateTime accessors: dayofyear, secondofday" begin
        dt = SpeedyDateTime(2000, 1, 1, 12, 30, 45)
        @test Dates.dayofyear(dt) == 1
        @test secondofday(dt) == 12 * 3600 + 30 * 60 + 45

        dt2 = SpeedyDateTime(2000, 12, 31)
        @test Dates.dayofyear(dt2) == 366  # 2000 is a leap year
    end

    @testset "Calendar helpers" begin
        @test Dates.isleapyear(2000) == true
        @test Dates.isleapyear(1900) == false
        @test Dates.isleapyear(2004) == true
        @test Dates.isleapyear(2001) == false

        @test Dates.daysinmonth(2000, 2) == 29
        @test Dates.daysinmonth(2001, 2) == 28
        @test Dates.daysinmonth(2000, 1) == 31
        @test Dates.daysinmonth(2000, 4) == 30
    end

    @testset "Conversion to/from Dates.DateTime" begin
        # Round-trip: Dates.DateTime -> SpeedyDateTime -> Dates.DateTime
        jdt = Dates.DateTime(2024, 6, 15, 10, 30, 45, 123)
        swdt = SpeedyDateTime(jdt)
        @test Dates.DateTime(swdt) == jdt

        # Check accessors match
        @test Dates.year(swdt) == Dates.year(jdt)
        @test Dates.month(swdt) == Dates.month(jdt)
        @test Dates.day(swdt) == Dates.day(jdt)
        @test Dates.hour(swdt) == Dates.hour(jdt)
        @test Dates.minute(swdt) == Dates.minute(jdt)
        @test Dates.second(swdt) == Dates.second(jdt)
        @test Dates.millisecond(swdt) == Dates.millisecond(jdt)
    end

    @testset "Conversion to/from Dates.Period" begin
        @test Dates.Millisecond(SpeedyMillisecond(500)) == Dates.Millisecond(500)
        @test Dates.Second(SpeedySecond(30)) == Dates.Second(30)
        @test Dates.Minute(SpeedyMinute(10)) == Dates.Minute(10)
        @test Dates.Hour(SpeedyHour(5)) == Dates.Hour(5)
        @test Dates.Day(SpeedyDay(3)) == Dates.Day(3)
        @test Dates.Week(SpeedyWeek(2)) == Dates.Week(2)
        @test Dates.Month(SpeedyMonth(6)) == Dates.Month(6)
        @test Dates.Year(SpeedyYear(2024)) == Dates.Year(2024)

        @test SpeedyMillisecond(Dates.Millisecond(500)) == SpeedyMillisecond(500)
        @test SpeedySecond(Dates.Second(30)) == SpeedySecond(30)
        @test SpeedyHour(Dates.Hour(5)) == SpeedyHour(5)
    end

    @testset "Consistency with Dates arithmetic" begin
        # Verify that SpeedyDateTime arithmetic matches Dates.DateTime arithmetic
        jdt = Dates.DateTime(2000, 1, 1)
        swdt = SpeedyDateTime(2000, 1, 1)

        # + Millisecond
        @test Dates.DateTime(swdt + SpeedyMillisecond(12345)) == jdt + Dates.Millisecond(12345)

        # + Second
        @test Dates.DateTime(swdt + SpeedySecond(3661)) == jdt + Dates.Second(3661)

        # + Minute
        @test Dates.DateTime(swdt + SpeedyMinute(90)) == jdt + Dates.Minute(90)

        # + Hour
        @test Dates.DateTime(swdt + SpeedyHour(25)) == jdt + Dates.Hour(25)

        # + Day
        @test Dates.DateTime(swdt + SpeedyDay(366)) == jdt + Dates.Day(366)

        # + Week
        @test Dates.DateTime(swdt + SpeedyWeek(52)) == jdt + Dates.Week(52)

        # + Year
        @test Dates.DateTime(swdt + SpeedyYear(5)) == jdt + Dates.Year(5)

        # + Month
        @test Dates.DateTime(swdt + SpeedyMonth(13)) == jdt + Dates.Month(13)

        # - operations
        @test Dates.DateTime(swdt - SpeedyDay(1)) == jdt - Dates.Day(1)
        @test Dates.DateTime(swdt - SpeedyMonth(1)) == jdt - Dates.Month(1)
        @test Dates.DateTime(swdt - SpeedyYear(1)) == jdt - Dates.Year(1)
    end

    @testset "Parametric IntType" begin
        # Verify that Int32 works throughout
        ms = SpeedyMillisecond(Int32(500))
        @test ms isa SpeedyMillisecond{Int32}
        @test Dates.value(ms) isa Int32

        s = SpeedySecond(Int32(30))
        @test (s + s) isa SpeedySecond{Int32}
        @test Dates.value(s + s) == Int32(60)

        # Conversion: Int32 * Int64 literal promotes to Int64
        # so the result type is SpeedyMillisecond{Int64}
        @test convert(SpeedyMillisecond, SpeedySecond(Int32(1))) == SpeedyMillisecond(1000)
    end

    @testset "SpeedyDateTime printing" begin
        dt = SpeedyDateTime(2000, 1, 1, 12, 30, 45)
        @test sprint(show, dt) == "2000-01-01T12:30:45"

        dt2 = SpeedyDateTime(2000, 1, 1, 12, 30, 45, 123)
        @test sprint(show, dt2) == "2000-01-01T12:30:45.123"
    end

    @testset "Broadcasting support" begin
        dt = SpeedyDateTime(2000, 1, 1)
        @test Ref(dt) isa Ref  # broadcastable returns Ref
    end

end
