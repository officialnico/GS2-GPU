!> A program that repeatedly calls an iglo loop for benchmarking

! Copyright Joseph Parker 2019 (joseph.parker@stfc.ac.uk)

! This file is part of gs2.

! Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files
! (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify,
! merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:

! The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
! LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
! CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

program test_indexing
  use unit_tests
  use benchmarks, only: benchmark_identifier
  use file_utils, only: append_output_file, close_output_file
  use gs2_init, only: init, init_level_list
  use gs2_layouts, only: g_lo, xxf_lo, yxf_lo, le_lo, layout
  use gs2_layouts, only: is_idx, it_idx, ik_idx, il_idx, ie_idx, ig_idx, isign_idx, idx
  use gs2_main, only: gs2_program_state_type, initialize_gs2, finalize_gs2
  use job_manage, only: time_message
  use mp, only: init_mp, finish_mp, proc0, nproc, mp_comm
  use run_parameters, only: nstep

  implicit none
  type(gs2_program_state_type) :: state
  real :: time_taken(2) = 0.0
  real :: time_init(2) = 0.0
  integer :: timing_unit

  ! General config
  ! (none)

  ! Set up depenencies
  call init_mp

  call announce_module_test('indexing')

  state%mp_comm_external = .true.
  state%mp_comm = mp_comm

  call initialize_gs2(state)

  ! initialize up to collisions terms level so all indexing is initialized
  call init(state%init, init_level_list%collisions)

  ! Begin test
  call announce_test('test index')
  call process_test(test_index(.true.), 'test index')

  call close_module_test('indexing')

  ! uninitialize
  call init(state%init, init_level_list%basic)
  call finalize_gs2(state)

  call close_module_test('test_indexing')

  call finish_mp

  contains
    function test_index(dummy)

      use unit_tests, only: announce_check, process_check, agrees_with
      use species, only: nspec
      use mp, only: proc0
      logical, intent(in) :: dummy
      logical :: test_index
      logical :: test_result   !! Overall test result, true if all pass
      logical :: check_result  !! Result of individual tests
      logical :: accelerated_x, accelerated_v
      character(2) :: istr

      test_result = .true.

      call announce_check('g layout')
      call check_glo(check_result)
      call process_check(test_result, check_result, 'g layout')

      call announce_check('xxf layout')
      call check_xxf(check_result)
      call process_check(test_result, check_result, 'xxf layout')

      call announce_check('yxf layout')
      call check_yxf(check_result)
      call process_check(test_result, check_result, 'yxf layout')

      call announce_check('le layout')
      call check_le(check_result)
      call process_check(test_result, check_result, 'le layout')

      test_index = test_result

    end function test_index

    subroutine check_glo(check_result)
      implicit none
      logical, intent(out) :: check_result
      integer :: i, j, k, l, it, ik, il, ie, is, iglo

      check_result = .true.

      do iglo=g_lo%llim_world, g_lo%ulim_world
        it=it_idx(g_lo,iglo)
        ik=ik_idx(g_lo,iglo)
        il=il_idx(g_lo,iglo)
        ie=ie_idx(g_lo,iglo)
        is=is_idx(g_lo,iglo)
        i=idx(g_lo,ik,it,il,ie,is)
        if(i.ne.iglo) then
          check_result = .false.
          exit
        end if
      end do

    end subroutine check_glo

    subroutine check_xxf(check_result)
      implicit none
      logical, intent(out) :: check_result
      integer :: i, ik, il, ie, is, ig, isgn, ixxf

      check_result = .true.

      do ixxf=xxf_lo%llim_world, xxf_lo%ulim_world
        ik=ik_idx(xxf_lo,ixxf)
        il=il_idx(xxf_lo,ixxf)
        ie=ie_idx(xxf_lo,ixxf)
        is=is_idx(xxf_lo,ixxf)
        ig=ig_idx(xxf_lo,ixxf)
        isgn=isign_idx(xxf_lo,ixxf)
        i=idx(xxf_lo,ig,isgn,ik,il,ie,is)
        if(i.ne.ixxf) then
          check_result = .false.
          exit
        end if
      end do

    end subroutine check_xxf

    subroutine check_yxf(check_result)
      implicit none
      logical, intent(out) :: check_result
      integer :: i, it, il, ie, is, ig, isgn, iyxf

      check_result = .true.

      do iyxf=yxf_lo%llim_world, yxf_lo%ulim_world
        it=it_idx(yxf_lo,iyxf)
        il=il_idx(yxf_lo,iyxf)
        ie=ie_idx(yxf_lo,iyxf)
        is=is_idx(yxf_lo,iyxf)
        ig=ig_idx(yxf_lo,iyxf)
        isgn=isign_idx(yxf_lo,iyxf)
        i=idx(yxf_lo,ig,isgn,it,il,ie,is)
        if(i.ne.iyxf) then
          check_result = .false.
          exit
        end if
      end do

    end subroutine check_yxf

    subroutine check_le(check_result)
      implicit none
      logical, intent(out) :: check_result
      integer :: i, it, ik, ig, is, ile

      check_result = .true.

      do ile=le_lo%llim_world, le_lo%ulim_world
        ig=ig_idx(le_lo,ile)
        it=it_idx(le_lo,ile)
        ik=ik_idx(le_lo,ile)
        is=is_idx(le_lo,ile)
        i=idx(le_lo,ig,ik,it,is)
        if(i.ne.ile) then
          check_result = .false.
          exit
        end if
      end do

    end subroutine check_le

end program test_indexing
