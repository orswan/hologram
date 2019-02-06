# hologram.jl
module hologram

using LinearAlgebra
using FFTW

#------------ SLM struct ----------------------------------------

struct SLM
	pixels
	pixelSize
	size
	voltLevels
	pixelx
	pixely
	centerx
	centery
	SLM(pixels,pixelSize,size,voltLevels) = new(pixels,pixelSize,size,voltLevels,
				range(pixelSize[1]/2,stop=size[1]-pixelSize[1]/2,length=pixels[1]),
				range(pixelSize[2]/2,stop=size[2]-pixelSize[2]/2,length=pixels[2]),
				size[1]/2,size[2]/2)
end

SLM(pixels,pixelSize,voltLevels) = SLM(pixels,pixelSize,pixels.*pixelSize,voltLevels)

Stest = SLM((10,8),(2.4,3.4),range(0,stop=3.4pi,length=102))
Stest2 = SLM((10,8),(2.4,3.4),range(0,stop=2pi,length=52))
S1920 = SLM((1920,1152),(9.2,9.2),range(0,stop=2pi,length=4096))
S512 = SLM((512,512),(15,15),range(0,stop=2pi,length=50))

#------------ Phase discretization ------------------------------

function discretizePhase(S::SLM,phase::Number)
	# Discretizes and mods phase onto available SLM levels of S.
	return S.voltLevels[ findmin(abs.( mod(phase,2pi) .- mod.(S.voltLevels,2pi) ))[2] ]
end

function discretizePhase(S::SLM,phase::Array{<:Number})
	# Discretizes an array of phases
	return [S.voltLevels[findmin(abs.(mod(p,2pi) .- mod.(S.voltLevels,2pi)))[2]] for p in phase]
end

function pixelizePhase(S::SLM,phase::Function;offset=[0,0])
	pixelx = S.pixelx .- S.centerx .+ offset[1]
	pixely = S.pixely .- S.centery .+ offset[2]
	return [phase(i,j) for i in pixelx, j in pixely]
end

function discretizePhase(S::SLM,phase::Function;offset=[0,0])
	return discretizePhase(S,pixelizePhase(S,phase))
end

#function plt(p)
#	# A useful function for visualizing phases
#	xs = vec([i for i=1:size(p)[1], j=1:size(p)[2]])
#	ys = vec([j for i=1:size(p)[1], j=1:size(p)[2]])
#	plot(xs,ys,vec(p),seriestype=:scatter)
#end

#---------------- Fourier transform -----------------------------

function subdivide(a::Array{<:Number,2},refinement)
	# Subdivides an array so that each lattice point expands to "refinement" number of lattice points.
	return kron(a,ones(refinement,refinement))
end

function Efield(S::SLM,phase::Union{Array{<:Number,2},Function},waist::Number; refinement=1)
	# Computes the electric field for a given phase profile and incident laser waist size.
		# Computes Fourier transform of phase with Gaussian envelope of given waist.
		# Refinement specifies how finely to subdivide the pixels. 
	if (typeof(phase)<:Array{<:Number,2}) && (size(phase) != (S.pixels...,))
		throw(ArgumentError,"Phase and number of pixels inconsistent")
	end
	discPhase = discretizePhase(S,phase)
	r2 = (S.pixelx .- S.centerx).^2 .+ (S.pixely' .- S.centery).^2		# Radius squared
	return exp.(-r2/waist^2 .+ im*discPhase)
end

function ft(S::SLM,phase::Union{Array{<:Number,2},Function},waist::Number;wrap=true)
	# Computes Fourier transform of phase with Gaussian envelope of given waist.
	if wrap
		out = abs.(fft(Efield(S,phase,waist))).^2
		sx,sy = size(out)
		return circshift(out, [floor(sx/2),floor(sy/2)])
	else
		return abs.(fft(Efield(S,phase,waist))).^2
	end
end

#------------------ Lensing with the SLM --------------------------

function hololens(S::SLM,f::Number)
	
	
end


end