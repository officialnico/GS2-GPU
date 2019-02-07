#!groovy

def build(){
	sh """
		module load ${params.modules}
	    	export GK_SYSTEM=${params.GK_SYSTEM}
	    	export MAKEFLAGS=-IMakefiles USE_NETCDF= USE_FFT= USE_NEW_DIAG=
	    	export TESTEXEC='mpirun -np 6'

	    	make depend
	    	make test_make
	    	make -j
    	"""
}

return this
