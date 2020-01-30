
g77omp = gfortran -fopenmp -fimplicit-none -fbackslash -O3


all: MD.IS-SPA.pol.sym.f
	$(g77omp) -o MD.IS-SPA.pol.sym.e MD.IS-SPA.pol.sym.f

md: MD.IS-SPA.pol.sym.f
	$(g77omp) -o MD.IS-SPA.pol.sym.e MD.IS-SPA.pol.sym.f

clean: 
	[ -f MD.IS-SPA.pol.sym.e ] && rm MD.IS-SPA.pol.sym.e || echo ''
