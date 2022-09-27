using JuliaFormatter

function format_code()
    i = 0
    while i < 10
        if format(dirname(@__DIR__))
            if i == 0
                return "Formatting unchanged"
            else
                return "Formatting corrected"
            end
        end
        i += 1
    end
    return "Formatting failed"
end

@info("Formatting...")
@info(format_code())
