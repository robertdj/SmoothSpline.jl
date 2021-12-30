@recipe function plot(sr::SplineRegression)
    @series begin
        all_x = range(sr.xmin, sr.xmax; length = 200)
        y = [SmoothSpline.predict(sr, x) for x in all_x]

        all_x, y
    end
end
