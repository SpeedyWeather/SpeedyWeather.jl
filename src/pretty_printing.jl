Base.show(io::IO, P::PrognosticVariables) = print(io,
    UnicodePlots.heatmap(gridded(P.temp[:,:,end])[:,end:-1:1]'.-273.15))