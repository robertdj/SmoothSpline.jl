# Example 2.2 on page 54 in "The NURBS Book"
function example_spline(t, i)
    y = 0.0

    if i == 0
        if 0 <= t < 1
            y = (1 - t)^2
        end
    elseif i == 1
        if 0 <= t < 1
            y = 2*t - 1.5*t^2
        elseif 1 <= t < 2
            y = 0.5 * (2 - t)^2
        end
    elseif i == 2
        if 0 <= t < 1
            y = 0.5*t^2
        elseif 1 <= t < 2
            y = -t^2 + 3*t - 1.5
        elseif 2 <= t < 3
            y = 0.5 * (3 - t)^2
        end
    elseif i == 3
        if 1 <= t < 2
            y = 0.5 * (t - 1)^2
        elseif 2 <= t < 3
            y = -t^2 + 5*t - 5.5
        elseif 3 <= t < 4
            y = 0.5 * (4 - t)^2
        end
    elseif i == 4
        if 2 <= t < 3
            y = 0.5 * (t - 2)^2
        elseif 3 <= t < 4
            y = -1.5*t^2 + 10*t - 16
        end
    elseif i == 5
        if 3 <= t < 4
            y = (t - 3)^2
        elseif 4 <= t < 5
            y = (t - 5)^2
        end
    elseif i == 6
        if 4 <= t < 5
            y = 2 * (t - 4) * (5 - t)
        end
    elseif i == 7
        if 4 <= t <= 5
            y = (t - 4)^2
        end
    end

    return y
end

