
##### poly fitting #####
using CurveFit

function curve_fitting!(turns::Int64, y::Vector{Float64}, n_pol::Int64) #y = records.ex...
    x = 1:1:turns
    fit = curve_fit(Polynomial, x, y, n_pol)
    y0 = fit.(x)

    plot(x, y0, linewidth=2, label = L"$y_{fitting}$")
    
end

curve_fitting!(10000, records.ex, 3)
plot!(xlabel="# Turns", title=L"Transverse \enspace x-emittance")



##### gaussian fitting #####

using LsqFit

# Nornalized Gaussian model
gaussian(x, p) = p[1] * exp.(-((x .- p[2]).^2) / (2 * p[3]^2))

#Data
y = @view(fftx[1:turns÷2])
y = real(y)
maximum(y)

x = @view(fftf[1:turns÷2])
#x = [x for x in f if x >= 0]

# Initial parameter guess: [amplitude, mean, standard deviation]
p0 = [maximum(y), 0.22, 0.01]

#just to check
#plot(x,y)

#fitting
fit = curve_fit(gaussian, x, y, p0)

# Extraction fit parameters
fitted_params = fit.param #amp, mean, std
#println("Fitted parameters: ", fitted_params)

# Generate points for the fitted curve
x_fit = range(0.2, 0.25, length=100)
y_fit = gaussian(x_fit, fitted_params)

# Plot the fitted Gaussian curve
plot(x_fit, y_fit, label=L"Gaussian Fit BB", linewidth=2)

ylabel!(L"Amplitude")
xlabel!(L"Tune")
