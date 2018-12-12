# Job manage tests

GS2 exits before nstep timesteps have been completed if:

1. a file `<job_name>.stop` is created; or
2. the job's run time (in seconds) exceeds the value of `avail_cpu_time` minus
   a 5 min = 300 sec grace period to allow gs2 to exit cleanly. 

These tests check that jobs are aborted by these two mechanisms. In `test_job_manage.in`,
`avail_cpu_time` is set to 302 seconds, so jobs should exit after 2 seconds.

* **`test_job_manage_avail_cpu_time`** checks that the job has stopped after two seconds,
  but before ten seconds.  Without limit, the job takes around 20 seconds.
* **`test_job_manage_stop_file`** creates a stop file, runs gs2, and checks
  that gs2 finishes _before_ two seconds have elapsed.
  
