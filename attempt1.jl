using Sundials
using PyPlot

function generateGrid(xmax, ymax, numPoints)
	x = linspace(0.0, xmax, numPoints)
	y = linspace(0.0, ymax, numPoints)
	#z = linspace(0, zmax, numPoints)
	grid = zeros(numPoints^2, 2)
	count = 1
	for i = 1:length(x)
		for j = 1:length(y)
			#for k = 1:length(z)
				#println("cords are $x[$i], $y[$j],  $z[$k]")
				
				grid[count,1] = x[i]
				grid[count,2] = y[j]
				#grid[count, 3] = z[k]
				count = count+1
			#end
		end
	end	

	return grid
end

function calculateGammaDot(u, v, deltaZ, deltaR, r)
	gammaDot = sqrt((v/deltaR)^2 +(v/r)^2 + (u/deltaZ)^2 + (u/deltaZ+u*deltaR)^2)
	return gammaDot
end

function calculateMu(gammaDot)
	#using Carreau model 
	#from http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1552100/

	mu0 = .056 #Pa s
	lambda = 3.313 # s
	n = .2568
	
	mu = mu0+2*mu0*(1+(lambda*gammaDot)^2)^(n-1)/2
	return mu
end

function momentum(t, w, wdot, p)
	#deltaX, deltaY, deltaZ, x, y
	#p is for passing in parameters
	deltaX = p[1]
	deltaY = p[2]
	deltaZ = p[3]
	x = p[4]
	y = p[5]

	v = w[1]
	u = w[2]

	P = 1 #need to actually deal with this later
	rho = 1.05 #g/cm^3, given by GENERALIZED APPROACH TO THE MODELING OF A R T E R I A L BLOOD FLOW

	deltaR = sqrt(deltaX^2 + deltaY^2)
	r = sqrt(x^2+y^2)
	gammaDot = calculateGammaDot(u, v, deltaZ, deltaR, r)
	mu = calculateMu(gammaDot)

	tauRR = -2*mu*v/deltaR
	tauZZ = -2*mu*u/deltaZ
	tauRZ = -1*mu*(u/deltaR + v/deltaZ)
	

	dvdt = -1*(u*v/deltaZ+v*v/deltaR)-1/rho*P/deltaR - 1/rho*(1/r*1/deltaR*(r*tauRR)+1/deltaZ*tauRZ)
	dudt = -1*(v*u/deltaR+u*u/deltaZ)-1/rho*P/deltaZ - 1/rho*(1/r*1/deltaR*(r*tauRZ)-1/deltaZ*tauZZ)
	wdot = [dvdt, dudt]
	return wdot
end
	
function main()
	#println("In main")
	numPoints = 100
	grid = generateGrid(1,1,numPoints)
	print(length(grid))
	deltat = .1;
	deltaX = grid[1,1]-grid[1,2]
	deltaY = grid[2,1]-grid[2,2]
	deltaZ = .1
	
	tFinal = 1;

	u0 = 1.0
	v0 = 1.0
	p0 = 1

	t0 = 0.0
	tsim = t0
	allData = Array[[ceil(tFinal/deltat)],[ceil(tFinal/deltat)]] 
	u = zeros(numPoints, numPoints)
	v = zeros(numPoints, numPoints)

	initials = [u0, v0]
	numRuns = 1

	while tsim < tFinal
		for j =1:numPoints-1
			for k= 1:numPoints-1
				xcord = grid[j,1]
				ycord = grid[k,2]
				println("At cordinates $xcord, $ycord")
				p = [deltaX, deltaY, deltaZ, xcord, ycord]
				wrappedMomentum(t,w,wdot) = momentum(t,w,wdot, p)
				tspan = [t0:deltat/10:t0+deltat]
				#println(typeof(tspan)) #they're floats, like they should be
				res = Sundials.cvode(wrappedMomentum, initials, tspan)
				#print("here")
				#vAll = float([ a[1] for a in res])
				#uAll = float([ a[2] for a in res])
				uAll = res[:,1]
				vAll = res[:,2]
				#get last elements (the ones at tf)
				u[j, k] = uAll[length(uAll)]
				#println("u is $u at $xcord, $ycord")
				v[j, k] = vAll[length(vAll)]
			end
		end
#	allData[numRuns, 1] = u
#	allData[numRuns, 2] = v		
	tsim = tsim+deltat
	numRuns = numRuns +1
	
	figure()
	pcolormesh(u)
	colorbar()
	title("z velocity")
	pcolormesh(v)
	colorbar()
	title("R velocity")
	end
end

main()
