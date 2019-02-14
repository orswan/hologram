# hologram.jl
module hologram

using LinearAlgebra
using FFTW

#------------ Typical usage -------------------------------------
"""
> import hologram
> using Plots
> plotly()
> H = hologram
> intensityProfile = abs.(H.ft(H.S512,( (x,y)-> (x^2 + y^2)/3e5 ),400;refinement=10)).^2
> plot(intensityProfile,st=:surface)
> intensityProfile = abs.(H.ft(H.S512,( (x,y)-> (x^2 + y^2)/3e4 + x/10 ),400;refinement=10)).^2;
> plot(intensityProfile,st=:surface)
"""
#------------ SLM struct ----------------------------------------

function midpointLattice(low::Number,high::Number,npoints::Integer)
	# Makes an array of npoints located symmetrically along interval [low,high], excluding the enpoints.
	return range(low,stop=high,length=(2*npoints+1))[2:2:end]
end

struct SLM
	pixels
	pixelSize
	size
	voltLevels
	pixelx
	pixely
	centerx
	centery
	cpixelx
	cpixely
	SLM(pixels,pixelSize,size,voltLevels) = new(pixels,pixelSize,size,voltLevels,
				midpointLattice(0,size[1],pixels[1]),
				midpointLattice(0,size[2],pixels[2]),
				size[1]/2,size[2]/2,
				midpointLattice(-size[1]/2,size[1]/2,pixels[1]),
				midpointLattice(-size[2]/2,size[2]/2,pixels[2])
				)
end

SLM(pixels,pixelSize,voltLevels) = SLM(pixels,pixelSize,pixels.*pixelSize,voltLevels)

function onPixel(S::SLM,x,y)
	# Returns true if the point (x,y) lies on a pixel.
#	println(latticize(x,S.cpixelx))
#	println(x)
#	println(S.pixelSize[1])
#	println(abs(latticize(x,S.cpixelx)-x))
	return (abs(latticize(x,S.cpixelx)-x)<=S.pixelSize[1]/2) && (abs(latticize(y,S.cpixely)-y)<=S.pixelSize[2]/2)
end

Stest = SLM((10,8),(2.4,3.4),range(0,stop=3.4pi,length=102))
Stest2 = SLM((10,8),(2.4,3.4),range(0,stop=2pi,length=52))
Stest3 = SLM((10,8),(2.4,2.4),(70.0,68.0),range(0,stop=2pi,length=52))
S1920 = SLM((1920,1152),(9.2,9.2),range(0,stop=2pi,length=4096))
S512 = SLM((512,512),(15,15),range(0,stop=2pi,length=50))
S512bad = SLM((512,512),(15,15),(14000,14000),range(0,stop=2pi,length=50))

#------------ Phase discretization ------------------------------

function latticize(x::Number,r::AbstractRange)
	# Finds nearest point of a range r to a given point x.
	if r[2]-r[1] < 0		# If given a backwards range, invert it before proceeding.
		R = r[end:-1:1]
	else
		R = r
	end
	if x < R[1]				# x lies below the range.
		return R[1]
	elseif x > R[end]		# x lies above the range.
		return R[end]
	else					# x lies within the range.
		idx = Int(round( (x-R[1])/(R[2]-R[1]) ))+1
		return R[idx]
	end
end
function latticize(x::Number,r::AbstractRange,modwhat::Number)
	# Finds the nearest point to x in r modulo modwhat.
	if r[2]-r[1] < 0			# If given a backwards range, invert it before proceeding.
		R = r[end:-1:1]
	else
		R = r
	end
	xs = range(x + ceil((R[1]-x)/modwhat) * modwhat, step=modwhat, stop=R[end])
	ls = [latticize(i,R) for i in xs]
	return ls[ findmin( abs.(ls-xs) )[2] ]
end

function discretizePhase(S::SLM,phase::Number)
	# Discretizes and mods phase onto available SLM levels of S.
	 return S.voltLevels[ findmin(abs.( mod(phase,2pi) .- mod.(S.voltLevels,2pi) ))[2] ]
#	return latticize(phase,S.voltLevels,2pi)		# This has about the same speed, but is faster without the 2pi.
end

function discretizePhase(S::SLM,phase::Array{<:Number})
	# Discretizes an array of phases
	return [S.voltLevels[findmin(abs.(mod(p,2pi) .- mod.(S.voltLevels,2pi)))[2]] for p in phase]
#	return [latticize(p,S.voltLevels,2pi) for p in phase]		# This has about the same speed, but is faster without the 2pi.
end

function pixelizePhase(S::SLM,phase::Function;offset=[0,0],refinement=1)
	if refinement==1
		pixelx = S.pixelx .- S.centerx .+ offset[1]
		pixely = S.pixely .- S.centery .+ offset[2]
		return [phase(i,j) for i in pixelx, j in pixely]
	elseif refinement>1
		xs = midpointLattice(0,S.size[1],S.pixels[1]*refinement) .- S.centerx .+ offset[1]
		ys = midpointLattice(0,S.size[2],S.pixels[2]*refinement) .- S.centery .+ offset[2]
		out = zeros(length(xs),length(ys))
		pixelx = S.pixelx .- S.centerx .+ offset[1]
		pixely = S.pixely .- S.centery .+ offset[2]
		phases = [phase(i,j) for i in pixelx, j in pixely]
		for i=1:refinement
			for j=1:refinement
				if onPixel(S,xs[i],ys[j])
					out[i:refinement:end,j:refinement:end] = phases
				end
			end
		end
		return out
	else
		throw(DomainError(refinement, "refinement must be a positive integer"))
	end
end

function discretizePhase(S::SLM,phase::Function;offset=[0,0],refinement=1)
	return discretizePhase(S,pixelizePhase(S,phase,offset=offset,refinement=refinement))
end

"""
function plt(p;step=1)
	# A useful function for visualizing phases
	xs = vec([i for i=1:size(p)[1], j=1:size(p)[2]])
	ys = vec([j for i=1:size(p)[1], j=1:size(p)[2]])
	plot(xs[1:step:end],ys[1:step:end],vec(p)[1:step:end],seriestype=:scatter)
end
"""

#---------------- Fourier transform -----------------------------

function subdivide(a::Array{<:Number,2},refinement)
	# Subdivides an array so that each lattice point expands to "refinement" number of lattice points.
	return kron(a,ones(refinement,refinement))
end
function subdivide(a::Array{<:Number,1},refinement)
	# Subdivides vector in analogy to the above array method.
	return kron(a,ones(refinement))
end

function Efield(S::SLM,phase::Union{Array{<:Number,2},Function},waist::Number; refinement::Int=1)
	# Computes the electric field for a given phase profile and incident laser waist size.
		# Computes Fourier transform of phase with Gaussian envelope of given waist.
		# Refinement specifies how finely to subdivide the pixels. 
	if (typeof(phase)<:Array{<:Number,2}) && (size(phase) != (S.pixels...,))
		throw(ArgumentError,"Phase and number of pixels inconsistent")
	end
	#discPhase = subdivide(discretizePhase(S,phase),refinement)
	discPhase = discretizePhase(S,phase,refinement=refinement)
	xs = midpointLattice(0,S.size[1],S.pixels[1]*refinement)
	ys = midpointLattice(0,S.size[2],S.pixels[2]*refinement)
	r2 = (xs .- S.centerx).^2 .+ (ys' .- S.centery).^2		# Radius squared
	return exp.(-r2/waist^2 .+ im*discPhase)
end

function ft(S::SLM,phase::Union{Array{<:Number,2},Function},waist::Number;wrap=true,refinement::Int=1)
	# Computes Fourier transform of phase with Gaussian envelope of given waist.
	if wrap
		out = abs.(fft(Efield(S,phase,waist,refinement=refinement))).^2
		sx,sy = size(out)
		return circshift(out, [floor(sx/2),floor(sy/2)])
	else
		return abs.(fft(Efield(S,phase,waist,refinement=refinement))).^2
	end
end

#------------------ Lensing with the SLM --------------------------

function hololens(S::SLM,f::Number)
	
	
end



#------------------ Test code -------------------------------------

"""
import hologram
H = hologram
p(x,y) = (x.^2 + y.^2)/r0^2 + x/x0
x0 = 1000
r0 = 10
S = H.SLM((10,10),(4,4),(50,50),range(0,stop=2pi,length=50))

plot(H.pixelizePhase(S,p,refinement=10),st=:surface)
plot(H.discretizePhase(S,p,refinement=10),st=:surface)
plot(H.ft(S,p,20,refinement=10),st=:surface,c=:blues)
x0 = 3
plot(H.ft(S,p,20,refinement=10),st=:surface,c=:blues)

"""



end
