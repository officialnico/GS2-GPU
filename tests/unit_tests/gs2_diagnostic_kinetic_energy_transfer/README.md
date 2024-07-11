Folder with testing scripts and results for the gs2 diagnostics which calculate the kinetic energy transfer to zonal flows.

## Idea
Run a cheap nonlinear test, the input file here is provided as "test_kinetic_energy_transfer.in". 
This file is similar to the nonlinear CBC test run. Then run the python script "" for the output and test
it against the output from the GS2 diagnostic with the file "".

## Updating golden answers
Assuming one has run the test case and you wish to update the golden answers then running `python compute_kinetic_energy_transfer_golden_answers.py --file test_kinetic_energy_transfer.out.nc --it -2 --write` should do this.
