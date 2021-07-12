# Helper file to generate the boilerplate files included
# in gs2 init.


default_module_routine_string = """
subroutine {name}_subroutine(going_up)
  use unit_tests, only: debug_message
  use {module_name}, only: {init_routine}, {finish_routine}
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- {name} ({level})')
    call {init_routine}
    current%level = init_level_list%{name}
  else
    call debug_message(1, 'gs2_init:init down reached init level -- {name} ({level})')
    call {finish_routine}
    current%level = init_level_list%{name} - 1
  end if
end subroutine {name}_subroutine
"""

more_general_module_routine_string = """
subroutine {name}_subroutine(going_up)
  use unit_tests, only: debug_message
{imports}
  logical, intent(in) :: going_up
  if (going_up) then
    call debug_message(1, 'gs2_init:init up reached init level -- {name} ({level})')
{up_action}
    current%level = init_level_list%{name}
  else
    call debug_message(1, 'gs2_init:init down reached init level -- {name} ({level})')
{down_action}
    current%level = init_level_list%{name} - 1
  end if
end subroutine {name}_subroutine
"""

class module_with_dependencies():
    def __init__(self,
                 name, # The name of the level/module
                 module_name = None, # The name of the module
                 dependencies = [], # A list of any modules this module depends on, identified by name
                 init_routine = None, # The name of the routine which initialises this module
                 finish_routine = None, # The name of the routine which finishes this module
                 visit_routine = None, # Can be used to hard-code the generated method
                 visit_args = {}
    ):

        self.name = name
        if module_name is None:
            self.module_name = name
        else:
            self.module_name = module_name
        self.dependency_names = dependencies
        self.init_routine = init_routine
        self.finish_routine = finish_routine
        self.visit_routine = visit_routine
        self.visit_args = visit_args
        self.level = 0

        if self.init_routine is None:
            self.init_routine = "init_{name}".format(name = self.name)
        if self.finish_routine is None:
            self.finish_routine = "finish_{name}".format(name = self.name)
        if self.visit_routine is None:
            self.visit_routine = default_module_routine_string

    def generate_method(self):
        return self.visit_routine.format(
            name = self.name,
            module_name = self.module_name,
            init_routine = self.init_routine,
            finish_routine = self.finish_routine,
            level = self.level,
            **self.visit_args
            ).strip()

if __name__ == "__main__":
    # First construct the modules / levels
    levels = [
        module_with_dependencies("full",
                                 dependencies = [
                                     "set_initial_values", "normalisations"],
                                 visit_routine = more_general_module_routine_string,
                                 visit_args = {'imports': "",
                                               'up_action': "",
                                               'down_action':""}
        ),
        module_with_dependencies("fields_level_2",
                                 module_name = "fields",
                                dependencies = [
                                    "collisions", "antenna",
                                    "dist_fn_level_3", "fields_level_1"]),
        module_with_dependencies("dist_fn_level_3",
                                 module_name = "dist_fn",
                                 dependencies = [
                                     "dist_fn_level_2", "hyper", "override_timestep"]),
        module_with_dependencies("dist_fn_level_1",
                                 module_name = "dist_fn",
                                 dependencies = [
                                     "dist_fn_arrays"]),
        module_with_dependencies("fields_level_1",
                                 module_name = "fields",
                                 dependencies = [
                                     "kt_grids", "antenna", "gs2_layouts",
                                     "fields_parameters"]),
        module_with_dependencies("collisions",
                                 dependencies = [
                                     "species", "kt_grids", "gs2_layouts", "theta_grid",
                                     "le_grids", "dist_fn_layouts", "run_parameters",
                                     "dist_fn_level_3"]),
        module_with_dependencies("nonlinear_terms",
                                 dependencies = [
                                     "species", "kt_grids", "gs2_layouts", "theta_grid",
                                     "le_grids", "dist_fn_layouts",
                                     "override_optimisations"]),
        module_with_dependencies("gs2_layouts", dependencies = []),
        module_with_dependencies("set_initial_values",
                                 dependencies = [
                                     "fields_level_2", "init_g", "override_initial_values"],
                                 visit_routine = more_general_module_routine_string,
                                 visit_args = {'imports': "",
                                               'up_action': "    call set_initial_field_and_dist_fn_values(current)",
                                               'down_action':""}
        ),
        module_with_dependencies("le_grids",
                                 dependencies = [
                                     "species", "kt_grids", "gs2_layouts", "theta_grid",
                                     "override_profiles"],
                                 visit_routine = more_general_module_routine_string,
                                 visit_args = {'imports':"  use le_grids, only: init_le_grids, finish_le_grids\n  logical :: dummy1, dummy2",
                                               'up_action': "    call init_le_grids(dummy1, dummy2)",
                                               'down_action':"    call finish_le_grids"}
        ),
        module_with_dependencies("antenna",
                                 dependencies = [
                                     "species", "le_grids", "run_parameters",
                                     "override_profiles"]),
        module_with_dependencies("theta_grid",
                                 dependencies = [
                                     "theta_grid_params", "override_miller_geometry"]),
        module_with_dependencies("normalisations", dependencies = []),
        module_with_dependencies("theta_grid_params", dependencies = []),
        module_with_dependencies("fields_parameters", module_name = "fields",
                                 dependencies = []),
        module_with_dependencies("kt_grids_parameters",
                                 module_name = "kt_grids",
                                 dependencies = [
                                     "theta_grid"]),
        module_with_dependencies("kt_grids",
                                 dependencies = [
                                     "theta_grid", "kt_grids_parameters",
                                     "override_kt_grids"]),
        module_with_dependencies("gs2_save", dependencies = []),
        module_with_dependencies("run_parameters",
                                 dependencies = [
                                     "kt_grids"]),
        module_with_dependencies("hyper",
                                 dependencies = [
                                     "kt_grids", "gs2_layouts"]),
        module_with_dependencies("init_g",
                                 dependencies = [
                                     "gs2_layouts"]),
        module_with_dependencies("species",
                                 dependencies = [
                                     "kt_grids"]),
        module_with_dependencies("dist_fn_parameters",
                                 module_name = "dist_fn",
                                 dependencies = [
                                     "gs2_layouts", "species", "theta_grid", "kt_grids",
                                     "le_grids"]),
        module_with_dependencies("dist_fn_arrays",
                                 module_name = "dist_fn",
                                 dependencies = [
                                     "dist_fn_parameters", "run_parameters",
                                     "nonlinear_terms", "dist_fn_layouts", "le_grids"]),
        module_with_dependencies("dist_fn_layouts",
                                 dependencies = [
                                     "species", "le_grids", "kt_grids", "gs2_layouts",
                                     "theta_grid"],
                                 visit_routine = more_general_module_routine_string,
                                 visit_args = {'imports':"\n".join([
                                    "  use gs2_layouts, only: init_dist_fn_layouts, finish_dist_fn_layouts",
                                    "  use kt_grids, only: naky, ntheta0",
                                    "  use le_grids, only: nlambda, negrid",
                                    "  use species, only: nspec"]),
                                               'up_action': "    call init_dist_fn_layouts(naky, ntheta0, nlambda, negrid, nspec)",
                                               'down_action': "    call finish_dist_fn_layouts"}

        ),
        module_with_dependencies("dist_fn_level_2",
                                 module_name = "dist_fn",
                                 dependencies = [
                                     "dist_fn_level_1", "override_profiles", "le_grids"]),
        module_with_dependencies("override_kt_grids",
                                 dependencies = [
                                     "kt_grids_parameters"],
                                 visit_routine = more_general_module_routine_string,
                                 visit_args = {'imports': "  use kt_grids, only: ktso => set_overrides",
                                               'up_action': "    if (current%kt_ov%is_initialised()) call ktso(current%kt_ov)",
                                               'down_action':""}
        ),
        module_with_dependencies("override_optimisations",
                                 dependencies = [
                                     "gs2_layouts", "fields_parameters"],
                                 visit_routine = more_general_module_routine_string,
                                 visit_args = {'imports':"\n".join([
                                    "  use gs2_layouts, only: lso=>set_overrides",
                                    "  use fields, only: fso=>set_overrides",
                                    "  use dist_fn, only: dso=>set_overrides"]),
                                               'up_action':"\n".join([
                                    "    if (current%opt_ov%is_initialised()) call lso(current%opt_ov)",
                                    "    if (current%opt_ov%is_initialised()) call fso(current%opt_ov)",
                                    "    if (current%opt_ov%is_initialised()) call dso(current%opt_ov)"]),
                                               'down_action':""}
                                 ),
        module_with_dependencies("override_miller_geometry",
                                 dependencies = [
                                     "theta_grid_params"],
                                 visit_routine = more_general_module_routine_string,
                                 visit_args = {'imports': "  use theta_grid_params, only: tgpso => set_overrides",
                                               'up_action': "    if (current%mgeo_ov%is_initialised()) call tgpso(current%mgeo_ov)",
                                               'down_action':""}
        ),
        module_with_dependencies("override_profiles",
                                 dependencies = [
                                     "species"],
                                 visit_routine = more_general_module_routine_string,
                                 visit_args = {'imports':"\n".join([
                                    "  use dist_fn, only: dfso=>set_overrides",
                                    "  use species, only: sso=>set_overrides"]),
                                               'up_action':"\n".join([
                                    "    if (current%prof_ov%is_initialised()) call dfso(current%prof_ov)",
                                    "    if (current%prof_ov%is_initialised()) call sso(current%prof_ov)"]),
                                               'down_action':""}

        ),
        module_with_dependencies("override_timestep",
                                 dependencies = [
                                     "run_parameters"],
                                         visit_routine = more_general_module_routine_string,
                                 visit_args = {'imports': "  use run_parameters, only: rso => set_overrides",
                                               'up_action': "    if (current%tstep_ov%is_initialised()) call rso(current%tstep_ov)",
                                               'down_action':""}
        ),
        module_with_dependencies("override_initial_values",
                                 dependencies = [
                                     "fields_level_2", "init_g"],
                                 visit_routine = more_general_module_routine_string,
                                 visit_args = {'imports': "",
                                               'up_action': "",
                                               'down_action':""}
                                 )
        ]

    # Convert to a ordered dictionary to ease lookup (ordered to ensure
    # reproducibility)
    from collections import OrderedDict
    levels_dict = OrderedDict([(x.name, x) for x in levels])

    # Now we need to decide what order we have to initialise things in
    # by analysing the dependency chain. We don't really care about
    # being efficient here so can adopt a fairly brute force approach
    levels_satisfied = []
    logical_level = 0
    levels_remaining = list(levels_dict.keys())

    # General approach:
    # 1. Check all remaining levels to see which have all dependencies
    #    in level_assigned_modules.
    # 2. Add this subset to level_assigned_modules
    # 3. Repeat until no remaining modules
    while len(levels_remaining) > 0:
        logical_level += 1
        def all_dependencies_statisfied(module):
            return all([x in levels_satisfied for x in module.dependency_names])
        subset_satisfied = list(filter(all_dependencies_statisfied, [levels_dict[x] for x in levels_remaining]))

        levels_satisfied.extend([x.name for x in subset_satisfied])
        [levels_remaining.remove(x.name) for x in subset_satisfied]

        if logical_level > len(levels):
            print("Something appears broken in generate_gs2_init.py",
                  "\nWe've searched ",len(levels)," times but still have yet to",
                  "find all dependencies for ",levels_remaining)
            for case in levels_remaining:

               print("For ",case,"we have the following dependencies with found status:")
               print([(x in levels_satisfied,x) for x in levels_dict[case].dependency_names])
            raise Exception("Unable to find all dependencies")

    # Now we have a list of modules/levels in the order required (from left to right)
    # We can therefore update our modules with the relevant level. Take the opportunity
    # to construct an ordered dictionary in appropriate order
    ordered_levels = OrderedDict()
    for i, level in enumerate(levels_satisfied):
        # We start with levels at 2 as we always have the hard-coded
        # basic entry with level = 1
        levels_dict[level].level = i + 2
        ordered_levels[level] = levels_dict[level]

    # Now we can form the init level list entries
    init_level_list = "\n".join([
        "  integer :: {name} = {level}".format(
            name = x.name,
            level = x.level,
        ) for x in ordered_levels.values()])

    # Now we can form the up/down calls
    up_calls = "\n".join([
        "    if (up() .and. (current%level == init_level_list%{name} - 1)) call {name}_subroutine(.true.)".format(
            name = x.name) for x in ordered_levels.values()])

    down_calls = "\n".join([
        "    if (down() .and. (current%level == init_level_list%{name})) call {name}_subroutine(.false.)".format(
            name = x.name) for x in ordered_levels.values()][::-1])

    # Now the visit methods
    visit_methods = "\n\n".join([
        x.generate_method() for x in ordered_levels.values()])

    import os
    import sys
    src_dir = os.path.dirname(sys.argv[0])
    with open(os.path.join(src_dir,"gs2_init_level_list.inc"), "w") as output:
        output.write(init_level_list.strip())

    with open(os.path.join(src_dir,"gs2_init_up.inc"), "w") as output:
        output.write(up_calls.strip())

    with open(os.path.join(src_dir,"gs2_init_down.inc"), "w") as output:
        output.write(down_calls.strip())

    with open(os.path.join(src_dir,"gs2_init_subroutines.inc"), "w") as output:
        output.write(visit_methods.strip())
