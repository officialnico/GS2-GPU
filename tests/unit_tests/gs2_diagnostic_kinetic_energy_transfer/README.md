Folder with testing scripts and results for the gs2 diagnostics which calculate the kinetic energy transfer to zonal flows

## input file for testing

    use a cbc test input file.

## folders of test simulations:

    - nproc_equal_ntheta: ntheta = 12 -> 13, nproc = 13

    - nproc_larger_ntheta: ntheta = 12 -> 13, nproc = 16

    - nproc_smaller_ntheta_moduloNonZero: ntheta = 12 -> 13, nproc = 4

    - nproc_smaller_ntheta_moduloZero: ntheta = 14 -> 15, nproc = 5

## testing run

    gs2 is not on scratch on my viking so running job in scratch and copying the result here

## files in this folder

    - compute_Tv_theta_*.py: uses python post-run script to calculate transfer
    - wrapper_comupute_Tv_all: wrapper for the above python scripts
    - Tv_theta_*.npz: the output from the python files
    - test_kinetic_energy_transfer.in: the gs2 input file used for testing
    - job.slurm: the job script to run the test run
    - job_*.log: the output of the test file job
    - run_test.py: runs the test of comparing gs2 diagnostic result with the python script results
    - TEST_RESULT.txt: the output of the run_test.py script

## Note:

    - ran input file twice with same seed and the transfer results agreed exactly