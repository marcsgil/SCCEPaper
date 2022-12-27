fy(y,x,χ) = SA[-4χ*x[1], ( -exp(-x[2])*cos(0.5*y[1]) + exp(-2*x[2])*cos(y[1]) )/χ]
fx(y,x,χ) = SA[ (exp(-x[2])*sin(0.5*y[1]) - exp(-2*x[2])*sin(y[1]) )/(2χ), -χ*y[2]]
H(x,χ) = χ*x[1]^2+(1-exp(-x[2]))^2/(4χ)
H2(x,χ) = H(x,χ)^2 - (2exp(-2x[2])-exp(-x[2]))/4
coord_transformation((P,Q),χ) = SA[√(1-Q^2)*P/(2χ),-log(1-Q)]