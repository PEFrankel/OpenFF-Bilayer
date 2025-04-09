! 
!   UNICAMP - University of Campinas
!   School of Chemical Engeneering (FEQ)
!   Complex Systems Engineering Laboratory (LESC)
! 
!   Ph.D. Luiz Guilherme Lomonaco Germiniani 
!   Supervisor: Prof. PhD. Luis Fernando Mercier Franco
! 

PROGRAM ORDER

  IMPLICIT NONE

  CHARACTER (len=25) :: infile, ndxfile, s6file, ovtfile, hexfile, logfile
  CHARACTER (len=5)  :: atom, text
  INTEGER :: f, i, j, k, l
  INTEGER :: p, pivot, counter
  INTEGER :: molec, numb, step
  INTEGER :: n_frames, n_lipids, n_chains, atoms_per_lipid
  INTEGER :: first, last, bin, n_bins, n_elements, n_data, n_half
  INTEGER, ALLOCATABLE :: atoms_per_chain(:), indexes(:)
  REAL :: time
  REAL :: cpu_start, cpu_end, elapsed, average_time, run_time, finished  
  DOUBLE PRECISION :: e2e_vec_x, e2e_vec_y, e2e_vec_z
  DOUBLE PRECISION :: norm, angle, bin_factor, delta_theta, x_bin
  DOUBLE PRECISION :: box_x, box_y
  DOUBLE PRECISION :: sum_cos, sum_sin, dist_x, dist_y, dist_sq, rcut_sq, s6_i
  DOUBLE PRECISION :: Rax(30), Ray(30), Raz(30)
  DOUBLE PRECISION, ALLOCATABLE :: Rx(:), Ry(:), Rz(:)
  DOUBLE PRECISION, ALLOCATABLE :: hist(:), S6(:), sum_x(:), sum_y(:), frac(:)

! Program header
  WRITE(*,*) '#'
  WRITE(*,*) '#  UNICAMP - LESC'
  WRITE(*,*) '#  Programer: Luiz G. L. Germiniani'
  WRITE(*,*) '#'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!  USER INPUTS  !!!!!!!!!!!!!

! Program header
  WRITE(*,*) '#  This code calculate the hexagonal order parameter for the packing'
  WRITE(*,*) '#  of lipid chains and the tilt angle distribution.'
  WRITE(*,*) '#'

  READ(*,*) infile, n_frames, ndxfile, s6file, hexfile, ovtfile

! INPUT OVERWRITE
  !n_frames     = 186

! Standard files
  !infile   = 'traj_new.gro'
  !ndxfile  = 'ndxfile.txt'
  !s6file   = 's6.xvg'  
  !hexfile  = 'ovito.xyz'

! File inputs
  !WRITE(*,*) '>> Provide a trajectory file (.gro):'
  !READ(*,*) infile
  !WRITE(*,*) '>> How many frames are in the file?'
  !READ(*,*) n_frames
  !WRITE(*,*) '>> Provide a file with the indexes for atoms in the computed molecules:'
  !READ(*,*) ndxfile


  WRITE(*,*) '>  All set!'
  WRITE(*,*) '>  '
  WRITE(*,*) '> Starting...'

!!!!!!!!!!!!!  USER INPUTS  !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  CALL CPU_TIME(run_time)

! Opening files
  OPEN(2,file=ndxfile)
  OPEN(1,file=infile)
  WRITE(*,*) '> ', infile, 'opened!'
  WRITE(*,*) '> ', ndxfile, 'opened!'

  logfile = TRIM(infile) // '.log'
  OPEN(3,file=logfile)

! Retriving index data
  WRITE(*,*) '> Getting atoms selection info from the index file...'
  READ(2,'(a2,i5)') text, n_lipids
  READ(2,'(a2,i5)') text, n_chains
  ALLOCATE (atoms_per_chain(n_chains))
  DO i = 1, n_chains
    READ(2,'(a2,i5)') text, atoms_per_chain(i)
  END DO

! Allocating frame data storaging arrays
  atoms_per_lipid = 0
  DO i = 1, n_chains
    atoms_per_lipid = atoms_per_lipid + atoms_per_chain(i)
  END DO
  n_data = n_lipids*atoms_per_lipid

  ALLOCATE (indexes(n_data))
  ALLOCATE (Rx(n_data),Ry(n_data),Rz(n_data))

! Reading index file
  DO i = 1, n_data
    READ(2,'(i5,x,a5,x,i5,x,i5)') indexes(i), atom, molec, numb
  END DO

  CLOSE(2)
  WRITE(*,*) '> Done!'

! Data processing

  rcut_sq = 0.65d0*0.65d0
  ALLOCATE (S6(n_frames),frac(n_frames))
  n_elements = n_lipids*n_chains
  ALLOCATE (sum_x(n_elements),sum_y(n_elements))

  OPEN(7,file=hexfile)
  WRITE(7,'(i9)') n_elements/2
  WRITE(7,'(a1)') ' '
  OPEN(8,file=ovtfile)
  WRITE(8,'(i9)') n_elements/2
  WRITE(8,'(a1)') ' '

!!!!!!!!!!!!  THE MAIN LOOP  !!!!!!!!!!!!
  average_time = 0.0e0
  DO f = 1, n_frames
    
    WRITE(*,*) '>>> FRAME:', f
    CALL READ_FRAME(indexes, Rx, Ry, Rz, n_data, time, step, box_x, box_y)
        
    !WRITE(*,*) '>> Computing xy projections...'
    p = 1
    counter = 1
    DO l = 1, n_lipids
      DO k = 1, n_chains
        sum_x(counter) = 0.d0
        sum_y(counter) = 0.d0
        DO j = 1, atoms_per_chain(k) 
          sum_x(counter) = sum_x(counter) + Rx(p)
          sum_y(counter) = sum_y(counter) + Ry(p)
          p     = p+1
        END DO
        sum_x(counter) = sum_x(counter)/DBLE(atoms_per_chain(k)) - 0.5d0*box_x
        sum_y(counter) = sum_y(counter)/DBLE(atoms_per_chain(k)) - 0.5d0*box_y
        sum_x(counter) = sum_x(counter) - box_x*ANINT(sum_x(counter)/box_x)
        sum_y(counter) = sum_y(counter) - box_y*ANINT(sum_y(counter)/box_y)
        counter =  counter + 1
      END DO
    END DO

    !WRITE(*,*) '>> Computing S6...'    
    n_half = n_elements/2
    S6(f) = 0.d0
    frac(f)  = 0.d0 

    DO l = 1, n_half
      counter = 0
      sum_cos = 0.d0
      sum_sin = 0.d0
      DO k = 1, l-1
        dist_x = sum_x(k) - sum_x(l)
        dist_y = sum_y(k) - sum_y(l)
        dist_x = dist_x - box_x*ANINT(dist_x/box_x)
        dist_y = dist_y - box_y*ANINT(dist_y/box_y)
        dist_sq = dist_x*dist_x + dist_y*dist_y
        IF (dist_sq < rcut_sq) THEN
          counter = counter + 1
          norm = DSQRT(dist_sq)
          angle = DSIGN(DACOS(dist_x/norm),dist_y)
          sum_cos = sum_cos + DCOS(6.d0*angle) 
          sum_sin = sum_sin + DSIN(6.d0*angle) 
        END IF
      END DO
      DO k = l+1, n_half
        dist_x = sum_x(k) - sum_x(l)
        dist_y = sum_y(k) - sum_y(l)
        dist_x = dist_x - box_x*ANINT(dist_x/box_x)
        dist_y = dist_y - box_y*ANINT(dist_y/box_y)
        dist_sq = dist_x*dist_x + dist_y*dist_y
        IF (dist_sq < rcut_sq) THEN
          counter = counter + 1
          norm = DSQRT(dist_sq)
          angle = DSIGN(DACOS(dist_x/norm),dist_y)
          sum_cos = sum_cos + DCOS(6.d0*angle) 
          sum_sin = sum_sin + DSIN(6.d0*angle) 
        END IF
      END DO
      s6_i = DSQRT(sum_cos*sum_cos + sum_sin*sum_sin) / 6.d0
      WRITE(7,*) sum_x(l), sum_y(l), 0.d0, 0.1d0, 0.1d0, 0.0d0, s6_i, counter
      IF(f .EQ. n_frames) THEN
        WRITE(8,*) sum_x(l), sum_y(l), 0.d0, 0.2d0, 0.2d0, 0.0d0, s6_i, counter
      ENDIF
      IF(s6_i .GE. 0.72d0) THEN
        frac(f) = frac(f) + 1.0d0
      ENDIF
      S6(f) = S6(f) + s6_i
    END DO

    DO l = n_half+1, n_elements
      counter = 0
      sum_cos = 0.d0
      sum_sin = 0.d0
      DO k = n_half+1, l-1
        dist_x = sum_x(k) - sum_x(l)
        dist_y = sum_y(k) - sum_y(l)
        dist_x = dist_x - box_x*ANINT(dist_x/box_x)
        dist_y = dist_y - box_y*ANINT(dist_y/box_y)
        dist_sq = dist_x*dist_x + dist_y*dist_y
        IF (dist_sq < rcut_sq) THEN
          counter = counter + 1
          norm = DSQRT(dist_sq)
          angle = DSIGN(DACOS(dist_x/norm),dist_y)
          sum_cos = sum_cos + DCOS(6.d0*angle) 
          sum_sin = sum_sin + DSIN(6.d0*angle) 
        END IF
      END DO
      DO k = l+1, n_elements
        dist_x = sum_x(k) - sum_x(l)
        dist_y = sum_y(k) - sum_y(l)
        dist_x = dist_x - box_x*ANINT(dist_x/box_x)
        dist_y = dist_y - box_y*ANINT(dist_y/box_y)
        dist_sq = dist_x*dist_x + dist_y*dist_y
        IF (dist_sq < rcut_sq) THEN
          counter = counter + 1
          norm = DSQRT(dist_sq)
          angle = DSIGN(DACOS(dist_x/norm),dist_y)
          sum_cos = sum_cos + DCOS(6.d0*angle) 
          sum_sin = sum_sin + DSIN(6.d0*angle) 
        END IF
      END DO
      s6_i = DSQRT(sum_cos*sum_cos + sum_sin*sum_sin) / 6.d0
      WRITE(7,*) sum_x(l), sum_y(l), 0.d0, 0.1d0, 0.1d0, 0.0d0, s6_i, counter
      IF(f .EQ. n_frames) THEN
        WRITE(8,*) sum_x(l), sum_y(l), 0.d0, 0.2d0, 0.2d0, 0.0d0, s6_i, counter
      ENDIF
      IF(s6_i .GE. 0.72d0) THEN
        frac(f) = frac(f) + 1.0d0
      ENDIF
      S6(f) = S6(f) + s6_i
    END DO
    S6(f) = S6(f)/DBLE(n_elements)
    frac(f) = frac(f)/DBLE(n_elements)
    WRITE(*,*) "S6 =", S6(f)
    WRITE(*,*) "X6 =", frac(f)
    WRITE(*,*) '> Frame completed!' 

  END DO
!!!!!!!!!!!!  END OF LOOP  !!!!!!!!!!!!

  OPEN(16,file=s6file)
  WRITE(*,*) '> S6 = ', sum(S6(:))/dble(n_frames)
  DO f = 1, n_frames
    WRITE(16,*) f, S6(f), frac(f)
  END DO

  DEALLOCATE (S6)
  DEALLOCATE (sum_x,sum_y)

  CLOSE(16)
  CLOSE(1)
  CLOSE(3)
  CLOSE(7)
  CLOSE(8)

  WRITE(*,*) '> Files done!' 

  CALL CPU_TIME(finished)
  run_time = finished - run_time
  WRITE(*,*) '> Overall CPU time:', run_time
  WRITE(*,*) '> check end of file:', Rx(n_data), Ry(n_data), Rz(n_data) 


  WRITE(*,*) '> END OF PROGRAM'
END PROGRAM ORDER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!     END     !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Uses the ndx data to identify the data points with the desired atoms
! position information and store it on the R arrays. 
  SUBROUTINE READ_FRAME(ndx, Rx, Ry, Rz, rsize, t, frame, lx, ly)

    CHARACTER (len=75) :: sys_name
    CHARACTER (len=50) :: info
    CHARACTER (len=5) :: resname, atom
    INTEGER :: i, j, k
    INTEGER :: n_particles, rsize
    INTEGER :: molec, numb, frame
    INTEGER :: ndx(rsize)
    DOUBLE PRECISION :: Rx(rsize), Ry(rsize), Rz(rsize)
    DOUBLE PRECISION :: lx, ly, lz
    REAL :: x, y, z, vx, vy, vz
    REAL :: t

    !WRITE(*,*) '>> Reading frame data... '

   !Initializing frame
    READ(1,'(a75)') sys_name
    READ(1,*) n_particles

   !Getting frame data
    j = INDEX(sys_name,'t= ')
    k = INDEX(sys_name,'step= ')

    IF ( j .NE. 0 .AND. k .NE. 0) THEN
      info = sys_name((j+3):(k-1))
      READ(info,*) t
      info = sys_name((k+6):75)
      READ(info,*) frame
    END IF


   !Readind lipid molecules (DPPC)
    k = 1
    DO i = 1, n_particles
      READ(1,'(i5,a5,a5,i5,3f8.3,3f8.4)') molec, resname, atom, numb, x, y, z, vx, vy, vz
      IF ( numb .EQ. ndx(k) ) THEN
        Rx(k) = DBLE(x)
        Ry(k) = DBLE(y)
        Rz(k) = DBLE(z)
        k = k+1
        !WRITE(*,'(i5,a5,a5,i5,3f8.3,3f8.4)') molec, resname, atom, numb, x, y, z, vx, vy, vz
      END IF
    END DO

   !End of frame
    READ(1,'(3f10.5)') lx, ly, lz
 
    !WRITE(*,*) '>> Done!' 

  END SUBROUTINE READ_FRAME


