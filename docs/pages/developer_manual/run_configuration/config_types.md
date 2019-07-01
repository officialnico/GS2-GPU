title: Configuration types and namelist input

The runtime setup of the code is determined by a user supplied _input
file_. In practice this typically contains a number of Fortran
namelists, each containing a set of key value pairs (this may be a
subset of all keys in the namelist, but it is an error if this
contains any undeclared keys).

@todo There are some undocumented features of our input parsing. For
example we can use `!include` syntax to include other files within an
input file as it is parsed. We should document these features.

Unlike some codes we do not read in the input namelists in a single
place in the code, instead preferring to keep the namelists (and
associated configuration variable declarations) local to the
appropriate module. Whilst this has advantages for scoping etc. it
means that similar code to deal with reading the inputs can appear
throughout much of the code. There are several common steps associated
with reading in a particular namelist:

1. Set the default values of _all_ input variables.
2. On processor 0 check the namelist is present in the input file using [[input_unit_exist]] and if it does exist use Fortran's native `read` to read the namelist.
3. We must then `broadcast` these values to all processors.

Historically these three stages have not always been consistently
presented in the code and it has been difficult to determine the
default values in all cases. Recognising the repetition involved in
the above three stages for reading the many namelists across the code,
coupled with a desire to make the defaults and associated
documentation easily discoverable the handling of the namelists has
recently been reworked. In particular a derived type is now used to
represent the set of input variables associated with a particular
namelist. To ensure consistent handling across all namelists these
derived types extend a common base `abstract class`
[[abstract_config_type]]. A further motivation for this change is to
allow users and developers to programmatically change the
configuration without writing new input files -- one may simply
provide a configuration object instead.

## The configuration type

By adopting a consistent approach it is possible to localise almost
all unique information regarding a particular namelist in the
configuration type. As such it is possible to generate all the
"boilerplate" associated with the common tasks of reading and
broadcasting the configuration variables. The process for doing this will be discussed shortly, but first the configuration type structure will be defined with the help of an example.

    !> Used to represent the input configuration of species
    type, extends(abstract_config_type) :: species_config_type
       ! namelist : species_knobs
       ! indexed : false     
       !> Definition not yet provided. Used in calculation of
       !> slowing-down distribution when electrons are adiabatic
       real :: me = -1.0
       !> Number of kinetic species to evolve.
       integer :: nspec = 2
       !> Definition not yet provided. Used in calculation of
       !> slowing-down distribution when electrons are adiabatic
       real :: zi_fac = -1.0
     contains
       procedure, public :: read => read_species_config
       procedure, public :: write => write_species_config
       procedure, public :: reset => reset_species_config
       procedure, public :: broadcast => broadcast_species_config
       procedure, public, nopass :: get_default_name => get_default_name_species_config
       procedure, public, nopass :: get_default_requires_index => get_default_requires_index_species_config
    end type species_config_type

The above demonstrates a small configuration type used in the species
module to represent the `species_knobs` namelist. This namelist has
two real valued inputs, `me` and `zi_fac`, and one integer input,
`nspec`. Note the Ford styled documentation preceeding each
declaration, this will be automatically parsed when the documentation
website is generated providing consistent documentation of the
inputs. Further we provide the default values in the declaration, this
ensures we can accurately document them and ensures we don't
accidentally work with uninitialised values if we forget to provide
one in the input file and haven't set the default. The contents below
the `contains` statement determines which routines provide the
required functionality (defined by the `abstract_config_type` base
class) -- this code can be/is auto-generated, see later. In addition
to the variable declaration, documentation and initialisation there
are two further parts of the configuration type that are required:

1. The comment `! namelist : species_knobs` defines the name of the namelist that this configuration type will read from, here it is `species_knobs`.
2. The comment `! indexed : false` determines if this namelist type is
a so-called _indexed_ namelist. In this instance it is not and
`species_knobs` is a _standard_ namelist. An _indexed_ namelist is one
where we allow/expect multiple instances of the same namelist which
differ in the input file just by an appended `_N` (with N some integer
index starting at 1). For example, the `species_parameters` namelist
is _indexed_ and the input file is expected to potentially contain
`species_parameters_1`, `species_parameters_2` etc.

### Generated code

Once the variable declarations and these two metadata comments are given it is possible to automatically generate the rest of the code that is required to work with this configuration object. A python script `scripts/config_genaration.py` is provided to perform this generation. To demonstrate how this works an example of adding a new namelist/configuration will be given below.

Suppose we want to add a new namelist called `performance_knobs` which
has two logical variables, `tune_performance` and `write_performance`,
and one real one, `unbalanced_fraction`, to the `performance` module
. We would first add the following to `performance.f90` (in fact this
could be in any temporary file for now).

    !> Used to represent the input configuration of species
    type, extends(abstract_config_type) :: performance_config_type
       ! namelist : performance_knobs
       ! indexed : false
       !> Do we try to tune the performance
       logical :: tune_performance = .false.
       logical :: write_performance = .false.
       !> What is the unbalance we allow
       real :: unbalanced_fraction = 0.5
    end type performance_config_type       

Note we have not provided a `contains` section and we have forgotten
to document `write_perfomance`. Now we can run the code generation
tool as follows

    python scripts/config_generation.py /path/to/performance.f90

This will then generate the file `performance_placeholders.f90` where
you ran the script. This will contain several sections of code than
can be copy-paste into `performance.f90` in order to provide all the
boilerplate. The first entry in the file is a canonical form of the configuration type decaration, in this case it would read:

    !> FIXME: Add documentation
    type, extends(abstract_config_type) :: performance_config_type
       ! namelist : performance_knobs
       ! indexed : False
       !> Do we try to tune the performance
       logical :: tune_performance = .false.
       !> What is the unbalance we allow
       real :: unbalanced_fraction = 0.5
       !> FIXME: Add documentation
       logical :: write_performance = .false.
     contains
       procedure, public :: read => read_performance_config
       procedure, public :: write => write_performance_config
       procedure, public :: reset => reset_performance_config
       procedure, public :: broadcast => broadcast_performance_config
       procedure, public :: get_default_name => get_default_name_performance_config
       procedure, public :: get_default_requires_index => get_default_requires_index_performance_config
    end type performance_config_type

The main differences here are:

1. The `contains` section has been populated.
2. The variables are declared in alphabetical order.
3. Any undocumented variables (and the type) will get a `!> FIXME: Add documentation` comment.

Following the new type definition there are some "free lines of code"
followed by a block comment containing `! Following is for the
config_type` and this is followed by a number of complete subroutines
and functions. These subroutines and functions provide the main
boilerplate code including all of the type-bound routines. These can
be simply copy-paste into the `performance` module (the convention is
to put these at the bottom of the module as they are uninteresting
generated code) -- in the future we may consider putting this
generated code in a sub-module. Returning to the lines referred to as
"free lines", these consist of several chunks. The first line is
simply a declaration designed to be used to define the module level
instance of this configuration type. In this case it reads

    type(performance_config_type) :: performance_config

The remaining free lines are designed to be used with the module's
`read_parameters` routine. This is the routine that is typically
responsible for reading and coordinating the input parameters. To
allow for users to override (or circumvent) the input file this
routine should be able to take a configuration instance. In other words we should have

    subroutine read_parameters(performance_config_in)

The next free line simply gives the declaration for this argument

    type(performance_config_type), intent(in), optional :: performance_config_in

The remaining free lines should then typically be the first exectuable statements in the `read_parameters` routine. Here they read

    if (present(config_in)) performance_config = performance_config_in

    call performance_config%init(name = 'performance_knobs', requires_index = .false.)

    ! Copy out internal values into module level parameters
    tune_performance = performance_config%tune_performance
    unbalanced_fraction = performance_config%unbalanced_fraction
    write_performance = performance_config%write_performance

    exist = config_internal%exist


The first simply deals with the optional argument. The next is the
statement that actually does the work of reading and broadcasting the
inputs. The remaining lines are then responsible for copying the
configuration values out of the configuration type and into the
module/local scope variables. In the long run we may wish to use the
configuration type memebers directly but for now we use this copying
to avoid having to change the majority of the code.