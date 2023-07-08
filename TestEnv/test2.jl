using QuasiMonteCarlo, Distributions,Plots
lb = [-Inf,-Inf]
ub = [Inf,Inf]
n = 2^5
d = 2

s = QuasiMonteCarlo.sample(n,lb,ub,GridSample())
s = QuasiMonteCarlo.sample(n,lb,ub,Uniform())
s = QuasiMonteCarlo.sample(n,lb,ub,SobolSample())
s = QuasiMonteCarlo.sample(n,lb,ub,LatinHypercubeSample())
s = QuasiMonteCarlo.sample(n,lb,ub,LatticeRuleSample())
s = QuasiMonteCarlo.sample(n,lb,ub,HaltonSample())
##
s = QuasiMonteCarlo.sample(2^12,lb,ub,SobolSample())
scatter(p,s[1,:],s[2,:],label=false,)