using SpeedyWeatherInternals.ProgressMeters
import SpeedyWeatherInternals.ProgressMeters: render_line, durationstring, eta

# a progress meter printing into a buffer with a fixed terminal width
progress_buffer(n; width = 80, color = false, kwargs...) =
    IOContext(IOBuffer(), :displaysize => (24, width), :color => color)
output_string(io::IOContext) = String(take!(io.io))

@testset "ProgressMeters" begin
    @testset "elements" begin
        p = Progress(200; output = progress_buffer(200))
        p.counter = 50
        p.tcurrent = p.tinit + 10      # 10 s for 50 iterations, 30 s remaining
        @test print_element(Percentage(), p) == " 25%"
        @test print_element(ETA(), p) == "ETA: 0:00:30"
        @test print_element(ElapsedTime(), p) == "0:00:10"
        @test print_element(Counter(), p) == "50"
        @test print_element(Bar(length = 8), p) == "━━      "
        @test print_element(Bar(length = 4, fill = '#', empty = '.'), p) == "#..."
        @test print_element(Bar(), p) == "━━━━━      " * " "^9    # 20 if printed on its own
        @test print_element("abc", p) == "abc"
        @test eta(p, p.tcurrent) ≈ 30

        # not 100% before finished, full bar and total time once finished
        p.counter = 199
        @test print_element(Percentage(), p) == " 99%"
        p.finished = true
        @test print_element(Percentage(), p) == "100%"
        @test print_element(Bar(length = 3), p) == "━━━"
        @test print_element(ETA(), p) == "Time: 0:00:10"

        # no division by zero
        p = Progress(0; output = progress_buffer(0))
        @test print_element(Percentage(), p) == "  0%"
        @test print_element(ETA(), p) == "ETA: N/A"
    end

    @testset "durationstring" begin
        @test durationstring(0) == "0:00:00"
        @test durationstring(3725.9) == "1:02:05"
        @test durationstring(2 * 86400 + 3 * 3600 + 4 * 60 + 5) == "2 days, 3:04:05"
        @test durationstring(12.5 * 86400) == "12.50 days"
        @test durationstring(Inf) == "N/A"
        @test durationstring(NaN) == "N/A"
        @test durationstring(-1) == "N/A"
    end

    @testset "line and separators" begin
        p = Progress(100; elements = (Percentage(), Bar(length = 4), ETA(), Counter(), ""), output = progress_buffer(100))
        p.counter = 50
        p.tcurrent = p.tinit + 5
        # single spaces next to the bar, separator otherwise, empty elements left out
        @test render_line(p) == " 50% ━━   ETA: 0:00:05, 50"
        p.description = "Run:"
        p.separator = " | "
        @test render_line(p) == "Run:  50% ━━   ETA: 0:00:05 | 50"

        # a bar of length 0 is left out, no double separators
        p = Progress(100; elements = (Percentage(), Bar(length = 0), Counter()), output = progress_buffer(100))
        @test render_line(p) == "  0%, 0"

        # a bar without length fills the remaining width (minus one column)
        for width in (40, 80)
            p = Progress(100; elements = (Percentage(), Bar(), ETA()), output = progress_buffer(100; width))
            @test textwidth(render_line(p)) == width - 1
        end
        # and is left out if there's no space
        p = Progress(100; elements = (Percentage(), Bar(), Counter()), output = progress_buffer(100; width = 8))
        @test render_line(p) == "  0%, 0"
        # lines are cut to the terminal width minus one column so they don't wrap
        p = Progress(100; elements = (Percentage(), Counter(), "αβγ"), output = progress_buffer(100; width = 11))
        @test render_line(p) == "  0%, 0, α"
        p.output = progress_buffer(100; width = 5)
        @test render_line(p) == "  0%"

        # custom elements
        struct Custom <: AbstractProgressElement end
        ProgressMeters.print_element(::Custom, p::Progress) = "custom $(p.counter)"
        p = Progress(10; elements = (Counter(), Custom()), output = progress_buffer(10))
        @test render_line(p) == "0, custom 0"
    end

    @testset "several lines" begin
        elements = (Percentage(), ETA(), NewLine(indent = 2), Counter(), "abc", NewLine(), NewLine(), "last")
        p = Progress(100; elements, output = progress_buffer(100))
        p.counter = 50
        p.tcurrent = p.tinit + 5
        # indented by the description and `indent`, lines without elements left out
        @test ProgressMeters.render_lines(p) == [" 50%, ETA: 0:00:05", "  50, abc", "last"]
        p.description = "Run:"
        @test ProgressMeters.render_lines(p) == ["Run:  50%, ETA: 0:00:05", "       50, abc", "     last"]
        @test render_line(p) == join(ProgressMeters.render_lines(p), '\n')
        # a bar fills the width of its line only
        p = Progress(100; elements = (Percentage(), Bar(), NewLine(), Counter()), output = progress_buffer(100; width = 20))
        lines = ProgressMeters.render_lines(p)
        @test textwidth(lines[1]) == 19 && lines[2] == "0"

        # redraws move the cursor up to the first line and overwrite both lines
        io = progress_buffer(3)
        p = Progress(3; dt = 0, elements = (Counter(), NewLine(), "second"), output = io)
        next!(p)
        @test output_string(io) == "\r1\e[K\nsecond\e[K"
        next!(p)
        @test output_string(io) == "\r\e[A2\e[K\nsecond\e[K"
        next!(p)
        @test output_string(io) == "\r\e[A3\e[K\nsecond\e[K\n"     # final lines end with a newline
    end

    @testset "drawing" begin
        # dt = 0: drawn on every step, final line with newline
        io = progress_buffer(3)
        p = Progress(3; dt = 0, elements = (Percentage(), ETA()), output = io)
        next!(p)
        line = output_string(io)
        @test startswith(line, "\r 33%, ETA: ") && endswith(line, "\e[K")
        next!(p)
        next!(p)
        line = output_string(io)
        @test occursin(r"\r100%, Time: \d:\d\d:\d\d\e\[K\n$", line)     # final line after the 67% one
        @test p.finished
        next!(p)        # nothing after the final line
        finish!(p)
        @test isempty(output_string(io))

        # color only if the output supports it
        io = progress_buffer(1; color = true)
        p = Progress(1; dt = 0, color = :blue, output = io)
        p.printed = true
        finish!(p)
        @test occursin(Base.text_colors[:blue], output_string(io))
        io = progress_buffer(1; color = false)
        p = Progress(1; dt = 0, output = io)
        p.printed = true
        finish!(p)
        @test !occursin("\e[3", output_string(io))

        # dt = Inf: never drawn, so also no final line
        io = progress_buffer(10)
        p = Progress(10; dt = Inf, output = io)
        for _ in 1:10
            next!(p)
        end
        finish!(p)
        @test isempty(output_string(io))
        @test p.counter == 10 && p.finished && !p.printed

        # disabled: counts but never draws
        io = progress_buffer(10)
        p = Progress(10; dt = 0, enabled = false, output = io)
        for _ in 1:5
            next!(p)
        end
        next!(p; step = 5)
        finish!(p)
        @test isempty(output_string(io))
        @test p.counter == 10

        # finish! draws the final line if drawn before, also when n is not reached
        io = progress_buffer(10)
        p = Progress(10; dt = 0, elements = (Percentage(),), output = io)
        next!(p)
        finish!(p)
        @test endswith(output_string(io), "\r100%\e[K\n")
        @test p.counter == 10

        # steps of several iterations
        io = progress_buffer(10)
        p = Progress(10; dt = 0, elements = (Counter(),), output = io)
        next!(p; step = 4)
        @test endswith(output_string(io), "\r4\e[K")
    end

    @testset "next! without redraw does not allocate" begin
        p = Progress(10^9; dt = Inf, output = progress_buffer(10^9))
        next!(p)
        @test (@allocated next!(p)) == 0
        p = Progress(10^9; enabled = false, output = progress_buffer(10^9))
        next!(p)
        @test (@allocated next!(p)) == 0
    end
end
