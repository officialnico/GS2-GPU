header = """
!> Wrapper routines for reading variables from netcdf file.
"""

module_string = """
{header}
module simpledataio_read
  use iso_fortran_env, only: real32, real64
  use simpledataio_write, only: real_imaginary_dimension_name
  implicit none
  private
  public :: read_variable_with_offset, read_variable

{read_with_offset_interface}
{read_interface}
contains
{read_with_offset}
{read}
end module simpledataio_read
"""

# Given a data type and dimension count determine the unique identifier
def get_basic_id(data_type, dimensions):

    return data_type.replace("(kind=","_").replace(")","")+"_"+str(dimensions)

def generate_read_method(data_type, dimensions, with_offset = True):
    method_string = """
  subroutine read_variable{method_type}_{id}(sfile, variable_name, val)
#ifdef FORTRAN_NETCDF
    use netcdf, only: nf90_get_var, nf90_strerror
#endif
    use simpledataio, only: number_of_unlimited_dimensions, number_of_dimensions, netcdf_inputs
    use simpledataio, only: sdatio_file, sdatio_int_kind, set_count, set_start
  {declarations}
  {body}
  end subroutine read_variable{method_type}_{id}
"""
    declarations_string = """
    type(sdatio_file), intent(in) :: sfile
    character(*), intent(in) :: variable_name
    {data_type},{dimensionality} intent(out) :: val
    integer(sdatio_int_kind), dimension(:), allocatable :: starts, counts, offsets
    integer :: fileid, varid, status, n2
    integer :: n, ndims
    {storage_array_decl}
"""
    body_string = """
#ifdef FORTRAN_NETCDF
    {get_n2}
    allocate(starts(n2))
    allocate(counts(n2))
    allocate(offsets(n2))

    call netcdf_inputs(sfile, variable_name, fileid, varid, starts, counts, offsets)
    status =  nf90_get_var(fileid, varid+1, &
         {storage_variable}, start=int(starts), count=int(counts))

    {copy_out}

    if (.not. status .eq. 0) write (*,*) 'Error reading variable: ', &
         variable_name, ', ',  nf90_strerror(status), 'starts=', starts, 'counts=', counts

    deallocate(starts, offsets, counts)
#endif
"""

    # Start handling specific case

    if data_type.startswith("complex"):
        is_complex = True
    else:
        is_complex = False

    if dimensions == 0 and not is_complex:
        is_scalar = True
    else:
        is_scalar = False

    # Get enough information for the declaration string
    if dimensions == 0:
        dimensionality = ""
    else:
        dimensionality = " dimension("+(":, "*(dimensions-1))+":),"

    # Decide how to declare any required temporary storage
    if is_complex:
        storage_array_decl = "real(kind = kind(val)), {shape} :: temp_storage".format(
            shape = "dimension(2"+"".join([", size(val, {ind})".format(ind=i+1) for i in range(dimensions)])+")")
    elif is_scalar:
        storage_array_decl = "{data_type}, dimension(1) :: temp_storage".format(data_type = data_type)
    else:
        storage_array_decl = ""

    # Now produce full declarations_string
    the_declarations = declarations_string.format(
        data_type = data_type,
        dimensionality = dimensionality,
        storage_array_decl = storage_array_decl).rstrip()

    # Now work on the main body. First handle getting n2
    if is_complex:
        extra_dimensions = 1
    else:
        extra_dimensions = 0

    if with_offset:
        get_n2 = """
    call number_of_unlimited_dimensions(sfile, variable_name, n)
    n2 = n+{dimensions}+{extra_dimensions}
    ndims = number_of_dimensions(sfile, variable_name)
    if (ndims /= n2) then
      write (*,*) "WARNING: The variable you pass to read_variable must have the same size &
           & and shape (excluding unlimited dimensions and the complex dimension) &
           & as the variable in the output file, regardless of what the values of starts &
           & and counts are. If you want to use a pass a different array shape to val &
           & you need to use read_variable_no_offset, which will behave in the default &
           & way for netcdf. &
           & You are probably about to encounter a segmentation fault. "
    end if
""".format(dimensions = dimensions, extra_dimensions = extra_dimensions).strip()

    else:
        get_n2 = "    n2 = number_of_dimensions(sfile, variable_name)".strip()

    # Now work out what variable + indexing we're reading into.
    storage_variable = "{variable}{indexing}"
    if is_scalar or is_complex:
        variable = "temp_storage"
    else:
        variable = "val"

    if with_offset and not is_scalar:
        indexing = "({first_offset}{other_offsets})"
        if is_complex:
            first_offset = ":"
        else:
            first_offset = "offsets(1):"
        other_offsets = "".join([
            ", offsets({ind}):".format(ind = i+1) for i in range(1, dimensions+extra_dimensions) ])

        indexing = indexing.format(first_offset = first_offset, other_offsets = other_offsets)
    else:
        indexing = ""
    storage_variable = storage_variable.format(variable = variable, indexing = indexing)

    # Now work out what, if any, copying out we need to do
    if is_scalar:
        copy_out = "val = temp_storage(1)"
    elif is_complex:
        indexing = "".join([", :" for i in range(dimensions)])
        copy_out = "val = cmplx( temp_storage(1{indexing}), temp_storage(2{indexing}))".format(indexing = indexing)
    else:
        copy_out = ""

    # Now we can form the body string
    the_body = body_string.format(
        get_n2 = get_n2,
        storage_variable = storage_variable,
        copy_out = copy_out
        ).rstrip()

    # Now get final bits
    if with_offset:
        method_type = "_with_offset"
    else:
        method_type = ""

    return method_string.format(
        method_type = method_type,
        id = get_basic_id(data_type, dimensions),
        declarations = the_declarations,
        body = the_body).rstrip()

if __name__ == "__main__":
    from itertools import product
    import argparse
    parser = argparse.ArgumentParser(
        description='Generate simpledataio code for reading.')
    parser.add_argument("output")
    args = parser.parse_args()

    data_types = ['real(kind=real32)', 'real(kind=real64)', 'integer',
                  'character', 'complex(kind=real32)', 'complex(kind=real64)']
    dimensions = [0, 1, 2, 3, 4, 5, 6]

    main_cases = list(product(data_types, dimensions))

    read_with_offset = "\n".join(
        [generate_read_method(*case, with_offset = True) for case in main_cases]
    )

    read = "\n".join(
        [generate_read_method(*case, with_offset = False) for case in main_cases]
    )

    read_with_offset_interface = """
  interface read_variable_with_offset
{contents}
  end interface read_variable_with_offset
""".format(contents =
           "\n".join(["    module procedure read_variable_with_offset_{id}".format(id = get_basic_id(*case)) for case in main_cases])).rstrip()

    read_interface = """
  interface read_variable
{contents}
  end interface read_variable
""".format(contents =
           "\n".join(["    module procedure read_variable_{id}".format(id = get_basic_id(*case)) for case in main_cases])).rstrip()


    with open(args.output, 'w') as output:
        output.write(module_string.format(
            header = header.strip(),
            read_with_offset_interface = read_with_offset_interface,
            read_interface = read_interface,
            read_with_offset = read_with_offset,
            read = read))
