program laplsolv
Use omp_lib
!-----------------------------------------------------------------------
! Serial program for solving the heat conduction problem 
! on a square using the Jacobi method. 
! Written by Fredrik Berntsson (frber@math.liu.se) March 2003
! Modified by Berkant Savas (besav@math.liu.se) April 2006
!-----------------------------------------------------------------------
  integer, parameter                  :: n=1000, maxiter=1000
  double precision,parameter          :: tol=1.0E-3
  double precision,dimension(0:n+1,0:n+1) :: T
  double precision,dimension(n)       :: tmp1,tmp2,tmp_end
  double precision                    :: error,local_error,x
  double precision                    :: t1,t0
  integer                             :: i,j,k,local_start,local_end,nt,me,k_g
  character(len=20)                   :: str
  
  ! Set boundary conditions and initial values for the unknowns
  T=0.0D0
  T(0:n+1 , 0)     = 1.0D0
  T(0:n+1 , n+1)   = 1.0D0
  T(n+1   , 0:n+1) = 2.0D0
  

  ! Solve the linear system of equations using the Jacobi method
  t0= omp_get_wtime()

  !$omp parallel shared(error, T, k_g) private(k, nt, local_error, me, local_start, local_size, tmp1, tmp2, tmp_end)

  me = omp_get_thread_num()
  nt = omp_get_num_threads() 
  call calculate_local_problem_size(n, nt, me, local_start, local_size)
  !print *, "Me: ", me, "Nt: ", nt, "Local_start: ", local_start, "Local_size: ", local_size

  do k=1,maxiter
  !$omp single	
     error=0.0D0
     k_g=k;
  !$omp end single nowait
     local_error=0.0D0
 
     tmp1=T(1:n,local_start-1)
     tmp_end=T(1:n,(local_start+local_size))

 !$omp barrier

 !Intermediate columns has no boundary values used by other threads
     do j=local_start, (local_start+local_size-2)
        tmp2=T(1:n,j)
        T(1:n,j)=(T(0:n-1,j)+T(2:n+1,j)+T(1:n,j+1)+tmp1)/4.0D0
        local_error=max(local_error,maxval(abs(tmp2-T(1:n,j))))
        tmp1=tmp2
     end do
     
  !Last column is a special case since tmp_end needs to be used
     tmp2=T(1:n,local_start+local_size-1)
     T(1:n,local_start+local_size-1)=(T(0:n-1,local_start+local_size-1)+T(2:n+1,local_start+local_size-1)+tmp_end+tmp1)/4.0D0
     local_error=max(local_error,maxval(abs(tmp2-T(1:n,local_start+local_size-1))))
     
  !$omp critical
     error = max(local_error,error)
  !$omp end critical
 
  !$omp barrier
     if (error<tol) then
        exit
     end if
  !$omp barrier
  end do

  !$omp end parallel
  
  t1= omp_get_wtime()

  write(unit=*,fmt=*) 'Time:',t1-t0,'Number of Iterations:',k_g
  write(unit=*,fmt=*) 'Temperature of element T(1,1)  =',T(1,1)

  ! Uncomment the next part if you want to write the whole solution
  ! to a file. Useful for plotting. 
  
  !open(unit=7,action='write',file='result.dat',status='unknown')
  !write(unit=str,fmt='(a,i6,a)') '(',N,'F10.6)'
  !do i=0,n+1
  !   write (unit=7,fmt=str) T(i,0:n+1)  
  !end do
  !close(unit=7)
end program laplsolv

!calculates number of lines and size for each subproblem stored in local_start and local_size
subroutine calculate_local_problem_size(n_f, nt_f, me_f, local_start, local_size)
	integer n_f,nt_f, linesize, rest, me_f, linestart, local_start, local_size
	linesize=n_f/nt_f
	rest=mod(n_f,nt_f)
       ! print *, 'Rest: ', rest, 'Linesize: ' ,linesize
	if(me_f<rest) then
		linesize = linesize + 1
		linestart = me_f*linesize + 1
	else
		linestart = rest*(linesize+1)+(me_f-rest)*linesize + 1
	end if
		
	local_start=linestart
	local_size=linesize
end subroutine calculate_local_problem_size

