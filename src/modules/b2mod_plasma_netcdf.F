      module b2mod_plasma_netcdf
      use b2mod_types
      use b2us_plasma
     & , only : B2State
      implicit none

      contains

#ifndef NO_CDF
      subroutine write_b2fplasma_netcdf(filename, nCv, nFc, ns, tim,
     & itim, state)
      implicit none
#     include <netcdf.inc>
      character*(*) filename
      integer :: nCv, nFc, ns, itim
      real (kind=R8) :: tim
      type(B2State) :: state
      integer :: ncid, status
      integer :: dim_nCv, dim_nFc, dim_ns, dim_idir, dim_one
      integer :: d_one(1), d_cv(1), d_cv_ns(2), d_cv_idir_ns(3)
      integer :: d_fc_idir(2), d_fc_idir_ns(3)
      integer :: timeid, itimid
      integer :: dna0id, dpa0id, hce0id, hci0id, hcibid
      integer :: vla0id, vsa0id
      integer :: fhtid, fheid, fhiid, fhnid, fhjid, fhmid, fhpid
      integer :: itim1(1)
      real (kind=R8) :: tim1(1)
      external check_cdf_status

      write(*,*) 'Writing b2fplasma.nc'
      status = nf_create(filename, NF_CLOBBER, ncid)
      call check_cdf_status(status)

      status = nf_def_dim(ncid, 'nCv', nCv, dim_nCv)
      call check_cdf_status(status)
      status = nf_def_dim(ncid, 'nFc', nFc, dim_nFc)
      call check_cdf_status(status)
      status = nf_def_dim(ncid, 'ns', ns, dim_ns)
      call check_cdf_status(status)
      status = nf_def_dim(ncid, 'idir', 2, dim_idir)
      call check_cdf_status(status)
      status = nf_def_dim(ncid, 'one', 1, dim_one)
      call check_cdf_status(status)

      d_one(1) = dim_one
      d_cv(1) = dim_nCv
      d_cv_ns(1) = dim_nCv
      d_cv_ns(2) = dim_ns
      d_cv_idir_ns(1) = dim_nCv
      d_cv_idir_ns(2) = dim_idir
      d_cv_idir_ns(3) = dim_ns
      d_fc_idir(1) = dim_nFc
      d_fc_idir(2) = dim_idir
      d_fc_idir_ns(1) = dim_nFc
      d_fc_idir_ns(2) = dim_idir
      d_fc_idir_ns(3) = dim_ns

      call define_double(ncid, 'time', 1, d_one, timeid)
      call define_int(ncid, 'itim', 1, d_one, itimid)
      call define_double(ncid, 'dna0', 2, d_cv_ns, dna0id)
      call define_double(ncid, 'dpa0', 2, d_cv_ns, dpa0id)
      call define_double(ncid, 'hce0', 1, d_cv, hce0id)
      call define_double(ncid, 'hci0', 1, d_cv, hci0id)
      call define_double(ncid, 'hcib', 2, d_cv_ns, hcibid)
      call define_double(ncid, 'vla0', 3, d_cv_idir_ns, vla0id)
      call define_double(ncid, 'vsa0', 2, d_cv_ns, vsa0id)
      call define_double(ncid, 'fht', 2, d_fc_idir, fhtid)
      call define_double(ncid, 'fhe', 2, d_fc_idir, fheid)
      call define_double(ncid, 'fhi', 2, d_fc_idir, fhiid)
      call define_double(ncid, 'fhn', 2, d_fc_idir, fhnid)
      call define_double(ncid, 'fhj', 2, d_fc_idir, fhjid)
      call define_double(ncid, 'fhm', 3, d_fc_idir_ns, fhmid)
      call define_double(ncid, 'fhp', 3, d_fc_idir_ns, fhpid)

      status = nf_enddef(ncid)
      call check_cdf_status(status)

      tim1(1) = tim
      itim1(1) = itim
      status = nf_put_var_double(ncid, timeid, tim1)
      call check_cdf_status(status)
      status = nf_put_var_int(ncid, itimid, itim1)
      call check_cdf_status(status)

      status = nf_put_var_double(ncid, dna0id, state%co%dna0)
      call check_cdf_status(status)
      status = nf_put_var_double(ncid, dpa0id, state%co%dpa0)
      call check_cdf_status(status)
      status = nf_put_var_double(ncid, hce0id, state%co%hce0)
      call check_cdf_status(status)
      status = nf_put_var_double(ncid, hci0id, state%co%hci0)
      call check_cdf_status(status)
      status = nf_put_var_double(ncid, hcibid, state%co%hcib)
      call check_cdf_status(status)
      status = nf_put_var_double(ncid, vla0id, state%co%vla0)
      call check_cdf_status(status)
      status = nf_put_var_double(ncid, vsa0id, state%co%vsa0)
      call check_cdf_status(status)

      status = nf_put_var_double(ncid, fhtid, state%dv%fht)
      call check_cdf_status(status)
      status = nf_put_var_double(ncid, fheid, state%dv%fhe)
      call check_cdf_status(status)
      status = nf_put_var_double(ncid, fhiid, state%dv%fhi)
      call check_cdf_status(status)
      status = nf_put_var_double(ncid, fhnid, state%dv%fhn)
      call check_cdf_status(status)
      status = nf_put_var_double(ncid, fhjid, state%dv%fhj)
      call check_cdf_status(status)
      status = nf_put_var_double(ncid, fhmid, state%dv%fhm)
      call check_cdf_status(status)
      status = nf_put_var_double(ncid, fhpid, state%dv%fhp)
      call check_cdf_status(status)

      status = nf_close(ncid)
      call check_cdf_status(status)
      return
      end subroutine write_b2fplasma_netcdf

      subroutine define_double(ncid, name, ndim, dims, varid)
      implicit none
#     include <netcdf.inc>
      integer :: ncid, ndim, dims(*), varid, status
      character*(*) name
      external check_cdf_status

      status = nf_def_var(ncid, name, NF_DOUBLE, ndim, dims, varid)
      call check_cdf_status(status)
      return
      end subroutine define_double

      subroutine define_int(ncid, name, ndim, dims, varid)
      implicit none
#     include <netcdf.inc>
      integer :: ncid, ndim, dims(*), varid, status
      character*(*) name
      external check_cdf_status

      status = nf_def_var(ncid, name, NF_INT, ndim, dims, varid)
      call check_cdf_status(status)
      return
      end subroutine define_int
#else
      subroutine write_b2fplasma_netcdf(filename, nCv, nFc, ns, tim,
     & itim, state)
      implicit none
      character*(*) filename
      integer :: nCv, nFc, ns, itim
      real (kind=R8) :: tim
      type(B2State) :: state

      write(*,*) 'Warning: ', trim(filename),
     & ' requested but B2.5 was compiled without NetCDF'
      return
      end subroutine write_b2fplasma_netcdf
#endif

      end module b2mod_plasma_netcdf
