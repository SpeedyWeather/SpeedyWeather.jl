using SpeedyWeatherInternals.Utils.TracableDates
import SpeedyWeatherInternals.Utils.TracableDates: value, toms, secondofday, days, format,
    FixedPeriod, UTInstant, AbstractTime, AbstractDateTime
using Test

@testset "Custom Dates" begin

    @testset "Period construction" begin
        @test value(Millisecond(1)) == 1
        @test value(Second(2)) == 2
        @test value(Minute(3)) == 3
        @test value(Hour(4)) == 4
        @test value(Day(5)) == 5
        @test value(Week(6)) == 6
        @test value(Month(7)) == 7
        @test value(Year(8)) == 8

        # parametric type is inferred from input
        @test Second(Int32(5)) isa Second{Int32}
        @test Second(5) isa Second{Int64}
    end

    @testset "Period printing" begin
        @test sprint(print, Second(1)) == "1 second"
        @test sprint(print, Second(2)) == "2 seconds"
        @test sprint(print, Day(1)) == "1 day"
        @test sprint(print, Day(3)) == "3 days"
        @test sprint(print, Hour(0)) == "0 hours"
    end

    @testset "Period zero/one/isfinite" begin
        @test zero(Second) == Second(0)
        @test zero(Day(5)) == Day(0)
        @test one(Hour) == 1
        @test isfinite(Minute(10))
        @test iszero(Second(0))
        @test !iszero(Second(1))
    end

    @testset "Period unary operations" begin
        @test -Second(5) == Second(-5)
        @test +Second(5) == Second(5)
        @test abs(Second(-3)) == Second(3)
        @test sign(Second(-3)) == -1
        @test sign(Second(3)) == 1
        @test signbit(Second(-1)) == true
        @test signbit(Second(1)) == false
    end

    @testset "Period same-type arithmetic" begin
        @test Second(3) + Second(2) == Second(5)
        @test Second(5) - Second(2) == Second(3)
        @test Second(6) / Second(3) == 2.0
        @test Second(10) / 2 == Second(5)
        @test div(Second(7), Second(3)) == 2
        @test rem(Second(7), Second(3)) == Second(1)
        @test mod(Second(7), Second(3)) == Second(1)
        @test Second(3) * 2 == Second(6)
        @test 2 * Second(3) == Second(6)
        @test Hour(2) * 3 == Hour(6)
    end

    @testset "Period comparison" begin
        @test Second(1) == Second(1)
        @test Second(1) != Second(2)
        @test Second(1) < Second(2)
        @test Second(2) > Second(1)
        @test isless(Second(1), Second(2))
    end

    @testset "Period conversions (FixedPeriod)" begin
        # Millisecond <-> Second
        @test convert(Millisecond, Second(1)) == Millisecond(1000)
        @test convert(Second, Millisecond(1000)) == Second(1)
        @test_throws InexactError convert(Second, Millisecond(999))

        # Second <-> Minute
        @test convert(Second, Minute(1)) == Second(60)
        @test convert(Minute, Second(60)) == Minute(1)

        # Minute <-> Hour
        @test convert(Minute, Hour(1)) == Minute(60)
        @test convert(Hour, Minute(60)) == Hour(1)

        # Hour <-> Day
        @test convert(Hour, Day(1)) == Hour(24)
        @test convert(Day, Hour(24)) == Day(1)

        # Day <-> Week
        @test convert(Day, Week(1)) == Day(7)
        @test convert(Week, Day(7)) == Week(1)

        # Multi-step conversions
        @test convert(Millisecond, Hour(1)) == Millisecond(3600000)
        @test convert(Second, Day(1)) == Second(86400)
        @test convert(Millisecond, Day(1)) == Millisecond(86400000)
    end

    @testset "Period conversions (OtherPeriod)" begin
        @test convert(Month, Year(1)) == Month(12)
        @test convert(Year, Month(12)) == Year(1)
        @test_throws InexactError convert(Year, Month(5))
    end

    @testset "Period promotion and cross-type arithmetic" begin
        # FixedPeriod promotion
        @test Second(1) + Millisecond(500) == Millisecond(1500)
        @test Minute(1) + Second(30) == Second(90)
        @test Hour(1) - Minute(30) == Minute(30)
        @test Day(1) + Hour(12) == Hour(36)

        # Cross-type comparison
        @test Second(60) == Minute(1)
        @test Hour(1) == Minute(60)
        @test Day(1) == Hour(24)

        # OtherPeriod promotion
        @test Year(1) + Month(6) == Month(18)
    end

    @testset "toms" begin
        @test toms(Millisecond(1)) == 1
        @test toms(Second(1)) == 1000
        @test toms(Minute(1)) == 60000
        @test toms(Hour(1)) == 3600000
        @test toms(Day(1)) == 86400000
        @test toms(Week(1)) == 604800000
    end

    @testset "Period construct-from-Period" begin
        @test Millisecond(Second(1)) == Millisecond(1000)
        @test Second(Minute(1)) == Second(60)
        @test Minute(Hour(1)) == Minute(60)
        @test Hour(Day(1)) == Hour(24)
        @test Day(Week(1)) == Day(7)
        @test Month(Year(1)) == Month(12)
    end

    @testset "DateTime construction" begin
        dt = DateTime(2000, 1, 1)
        @test year(dt) == 2000
        @test month(dt) == 1
        @test day(dt) == 1
        @test hour(dt) == 0
        @test minute(dt) == 0
        @test second(dt) == 0
        @test millisecond(dt) == 0

        dt2 = DateTime(2024, 3, 15, 14, 30, 45, 123)
        @test year(dt2) == 2024
        @test month(dt2) == 3
        @test day(dt2) == 15
        @test hour(dt2) == 14
        @test minute(dt2) == 30
        @test second(dt2) == 45
        @test millisecond(dt2) == 123
    end

    @testset "DateTime from Period types" begin
        dt = DateTime(Year(2000), Month(6), Day(15))
        @test year(dt) == 2000
        @test month(dt) == 6
        @test day(dt) == 15
    end

    @testset "DateTime equality and comparison" begin
        dt1 = DateTime(2000, 1, 1)
        dt2 = DateTime(2000, 1, 1)
        dt3 = DateTime(2000, 1, 2)

        @test dt1 == dt2
        @test dt1 != dt3
        @test dt1 < dt3
        @test dt3 > dt1
        @test isless(dt1, dt3)
    end

    @testset "DateTime - DateTime" begin
        dt1 = DateTime(2000, 1, 1)
        dt2 = DateTime(2000, 1, 2)
        @test dt2 - dt1 == Millisecond(86400000)
        @test dt1 - dt1 == Millisecond(0)
    end

    @testset "DateTime +/- FixedPeriod" begin
        dt = DateTime(2000, 1, 1, 0, 0, 0)

        @test dt + Millisecond(1) == DateTime(2000, 1, 1, 0, 0, 0, 1)
        @test dt + Second(1) == DateTime(2000, 1, 1, 0, 0, 1)
        @test dt + Minute(1) == DateTime(2000, 1, 1, 0, 1, 0)
        @test dt + Hour(1) == DateTime(2000, 1, 1, 1, 0, 0)
        @test dt + Day(1) == DateTime(2000, 1, 2)
        @test dt + Week(1) == DateTime(2000, 1, 8)

        @test dt - Day(1) == DateTime(1999, 12, 31)
        @test dt - Hour(1) == DateTime(1999, 12, 31, 23, 0, 0)

        # commutativity
        @test Day(1) + dt == dt + Day(1)
    end

    @testset "DateTime +/- Year" begin
        dt = DateTime(2000, 2, 29)  # leap day
        @test dt + Year(1) == DateTime(2001, 2, 28)  # clamp to end of Feb
        @test dt - Year(1) == DateTime(1999, 2, 28)
        @test dt + Year(4) == DateTime(2004, 2, 29)  # still leap year
    end

    @testset "DateTime +/- Month" begin
        dt = DateTime(2000, 1, 31)
        @test dt + Month(1) == DateTime(2000, 2, 29)  # clamp to Feb 29 (leap year)
        @test dt + Month(2) == DateTime(2000, 3, 31)

        dt2 = DateTime(2001, 1, 31)
        @test dt2 + Month(1) == DateTime(2001, 2, 28)  # non-leap year

        # wrap around year boundary
        dt3 = DateTime(2000, 11, 15)
        @test dt3 + Month(3) == DateTime(2001, 2, 15)
        @test dt3 - Month(12) == DateTime(1999, 11, 15)
    end

    @testset "DateTime accessors: dayofyear, secondofday" begin
        dt = DateTime(2000, 1, 1, 12, 30, 45)
        @test dayofyear(dt) == 1
        @test secondofday(dt) == 12 * 3600 + 30 * 60 + 45

        dt2 = DateTime(2000, 12, 31)
        @test dayofyear(dt2) == 366  # 2000 is a leap year
    end

    @testset "Calendar helpers" begin
        @test isleapyear(2000) == true
        @test isleapyear(1900) == false
        @test isleapyear(2004) == true
        @test isleapyear(2001) == false

        @test daysinmonth(2000, 2) == 29
        @test daysinmonth(2001, 2) == 28
        @test daysinmonth(2000, 1) == 31
        @test daysinmonth(2000, 4) == 30
    end

    @testset "Parametric IntType" begin
        # Verify that Int32 works throughout
        ms = Millisecond(Int32(500))
        @test ms isa Millisecond{Int32}
        @test value(ms) isa Int32

        s = Second(Int32(30))
        @test (s + s) isa Second{Int32}
        @test value(s + s) == Int32(60)

        # Conversion: Int32 * Int64 literal promotes to Int64
        # so the result type is Millisecond{Int64}
        @test convert(Millisecond, Second(Int32(1))) == Millisecond(1000)
    end

    @testset "DateTime printing" begin
        dt = DateTime(2000, 1, 1, 12, 30, 45)
        @test sprint(show, dt) == "2000-01-01T12:30:45"

        dt2 = DateTime(2000, 1, 1, 12, 30, 45, 123)
        @test sprint(show, dt2) == "2000-01-01T12:30:45.123"
    end

    @testset "Broadcasting support" begin
        dt = DateTime(2000, 1, 1)
        @test Ref(dt) isa Ref  # broadcastable returns Ref
    end

    @testset "format helper" begin
        dt = DateTime(2024, 3, 15, 14, 30, 45)
        @test format(dt, "yyyy-mm-dd HH:MM:0.0") == "2024-03-15 14:30:0.0"
    end

    @testset "days and firstdayofmonth helpers" begin
        dt = DateTime(2000, 3, 15, 12, 0, 0)
        fdm = firstdayofmonth(dt)
        @test year(fdm) == 2000
        @test month(fdm) == 3
        @test day(fdm) == 1
        @test days(dt - fdm) == 14
    end

end
