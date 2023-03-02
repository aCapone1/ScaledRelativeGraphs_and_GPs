using Plots
using LinearAlgebra
using FFTW

# n is the number of harmonics of the Fourier transform to take.
# Around 30 seems good.
function sat_offset(n=30)
        θ = 11.11
        start = 0
        step = 0.1
        max = 4

        input_sin = sin.((LinRange(0, 2*pi/θ, 30)))
        function sat(x)
                if x > 1
                        return 1
                elseif x < -1
                        return -1
                end
                return x
        end

       amplitudes = 3 .^(start:step:max)
       offsets = -10.75:1:10.75

       z = []
       for A1 = amplitudes
               for A2 = amplitudes
                       for o1 = offsets
                               for o2 = offsets

                                       fdy = fft(sat.(o1 .+ A1.*input_sin) - sat.(o2 .+ A2.*input_sin))[1:n]
                                 
                                       fdu = fft((o1 .+ A1.*input_sin) - (o2 .+ A2.*input_sin))[1:n]
        
                                       ip = (fdu⋅fdy)/length(fdu)
                                       input_norm = sqrt((fdu⋅fdu/length(fdu)))
                                       output_norm = sqrt((fdy⋅fdy/length(fdy)))
                                       angle_arg = real(ip/input_norm/output_norm)
        
                                       angle = (abs(angle_arg) > 1 ? 0 : acos(angle_arg))
                                       gain = output_norm/input_norm
        
                                       append!(z, gain*exp(im*angle))
                                       append!(z, gain*exp(-im*angle))
                               end
                       end
               end
       end
       p = scatter(real.(z), imag.(z), markersize=2, markerstrokewidth=0, legend = :none, markercolor = :cornflowerblue, fontfamily="serif-roman",  framestyle = :origin, aspect_ratio = 1, xlabel = "Re(z)", ylabel = "Im(z)")
        #savefig(p, "DF_SRG.png")
        display(p)
end
