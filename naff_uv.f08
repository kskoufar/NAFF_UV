!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! The following NAFF_UV program is an enhanced revision of JACQUES LASKAR's original NAFF that was written in Fortran 77. 
! NAFF_UV is written in Fortran 2008 and is a heavily modified version of the Michael Ehrlichman source code found in https://github.com/MichaelEhrlichman/FortNAFF (Fortran 90)
!
! Exept of the main program the following code consists of the NEEDED_STRUCTURES, READ_STORE_DATA and NAFF_UV modules.
! There is also a dipendence on the KIND_ACCURACY and BRENT_MNBRAK_METHODS modules that are in separate files and from
! the FFTW libraries.
!
! About the main module NAFF_UV.
!
! This module implements the core of the NAFF algorithm for calculating the spectra of periodic data.
!
! This module is useful when you need an accurate frequency-amplitude decomposition from only a small number of samples.
!
! Precision of the determined freqiencies goes as 1/N^4, where N is the number of data points.
!
! The current version works with 1D complex (D(:) = D1(:) + i D2(:)) or 1D real signale.
!
!
! The steps of NAFF_UV are:
! 1) Estimate peak fr_n in frequency spectrum using "windowed" data (wdata) as input for an interpolated FFTW.
! 	 Differend windows can be used : gaussina and hann_p where p an integer
! 2) Refine fr_n (out_freq) estimate by using optimizer to maximize the -1*<e^{i fr_n}|wdata> (<e^{i fr_n}|wdata>=int[wdata*e^{-i fr_n}]) with the help of mnbral and brent algorithms.
! 3) Refine amp_n (camp) estimate by using the projection <e^{i fr_n}|wdata> (<e^{i fr_n}|wdata>=int[wdata*e^{-i fr_n}]).
! 4) Remove e_1 = amp_n*e^{i fr_n} component from the data.
! 5) Repeat step 1 to estimate the new frequency component.
! 6) Repeat steps 2-3 to refine the new frequency and amplitude component.
! 7) Use Gram-Schmidt to build a basis function by orthogonalizing the new frequency component relative to the existing frequency components.
!
! For more details on naff -> Frequency map analysis and quasiperiodic decompositions Jacques Laskar
!-

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

module NEEDED_STRUCTURES

  use KIND_ACCURACY

  implicit none

  type input_data

    real (kind=rac) :: real_part
    real (kind=rac) :: imaginary_part

  end type input_data


  type harm_freq_ampl

    integer (kind=iac) :: harmonic_num
    real (kind=rac) :: frequency
    real (kind=rac) :: amplitude
    real (kind=rac) :: real_ampl
    real (kind=rac) :: imaginary_ampl

  end type harm_freq_ampl


end module NEEDED_STRUCTURES

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

module READ_STORE_DATA

  use, intrinsic :: iso_c_binding
  use KIND_ACCURACY
  use NEEDED_STRUCTURES

  implicit none

  character (len=1), private :: comment_symbol='#'

contains

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine GetLines (unit, file_name, n_lines, n_active_lines)

    implicit none

    integer , intent(in) :: unit
    character (len=*), intent(in) :: file_name
    integer (kind=iac) , intent(out) :: n_lines, n_active_lines
    integer :: io
    character (len=1) :: comment_check

    open(unit=unit,file=file_name,status='old', action='read')

    n_lines = 0
    n_active_lines = 0
    do

      read(unit=unit, fmt=*, iostat=io) comment_check

      if (io > 0) then

        print*, "ERROR - subroutine GetLines !!!"
        print*, "Some error occured while reading the file with unit: ", unit, " !!!"
        rewind(unit)
        return

      elseif (io == 0 .and. comment_check /= comment_symbol) then

        n_active_lines = n_active_lines + 1 ! all the none commented lines

        n_lines = n_lines + 1

      elseif (io == 0 .and. comment_check == comment_symbol) then

        n_lines = n_lines + 1

      else

        !#!print*, "subroutine GetLines: - The file with unit: ",unit ," is read normally."
        !#!print*, "It has ", n_lines, "lines, ", n_lines-n_active_lines, "of them are coments and the rest", &
        !       n_active_lines, "are active."
        
        exit

      endif

    enddo

    rewind (unit)

    close (unit)

  end subroutine GetLines

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function GetData (unit, file_name, number_of_initial_line, f_p_number_active_lines, data_step, input_data_type) result(data)

    implicit none

    integer, intent(in) :: unit
    character (len=*), intent(in) :: file_name
    integer (kind=iac), intent(in) :: number_of_initial_line, f_p_number_active_lines, data_step
    character (len=*), intent(in) :: input_data_type
    type (input_data), dimension(f_p_number_active_lines) :: data
    integer (kind=iac) :: jj, ii, nn
    integer :: io
    character (len=1000) :: buffer
    character (len=:), allocatable :: line_data

    open (unit=unit, file=file_name, status='old', action='read')

    jj = 1
    do while (jj<number_of_initial_line)

      read (unit=unit, fmt="(A)", iostat=io) buffer
      line_data = trim(buffer)

      if (io > 0) then

        print*, "ERROR - subroutine GetData !!!"
        print*, "Some error occured while reading the file", file_name, " with unit: ", unit, " !!!"
        rewind(unit)
        return

      else if (io == 0) then

        if (line_data(1:1) == comment_symbol) then

          ! print*, "The line: ", ii, " is skipt since it is a comment!!!"
          cycle

        else

          jj=jj+1

        endif

      endif

    enddo


    ii = 0
    nn = 1
    do while (nn<=f_p_number_active_lines)

      read (unit=unit, fmt="(A)", iostat=io) buffer
      line_data = trim(buffer)


      if (io > 0) then

        print*, "ERROR - subroutine GetData !!!"
        print*, "Some error occured while reading the file", file_name, " with unit: ", unit, " !!!"
        rewind(unit)
        return

      else if (io == 0) then

        if (line_data(1:1) == comment_symbol) then

          ! print*, "The line: ", ii, " is skipt since it is a comment!!!"
          cycle

        elseif (modulo(ii,data_step) /= 0) then

          ii = ii + 1
          cycle

        else

          if (input_data_type=='c') then
            read (line_data, fmt=*) data(nn)%real_part, data(nn)%imaginary_part
          elseif (input_data_type=='r') then
            read (line_data, fmt=*) data(nn)%real_part
            data(nn)%imaginary_part=0.0_rac
            !cheack it befor use the follow (ex. amp*congj(amp) not true)
            !data(nn)%imaginary_part=transfer((huge(1_rac)), 1.0_rac)
          else
            print*, "ERROR - subroutine GetData !!!"
            print*, "The type of the input data must be complex <c> or real <r> !!!"
            rewind(unit)
            return
          endif
          ii = ii+1
          nn = nn+1

        endif

      elseif (io < 0) then

        !#!print*, "function GetData: - The file ", file_name, " with unit: ",unit ," is read normally."

      endif

    enddo

    rewind(unit)
    close (unit)

  end function GetData

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine StoreFreqAmp (unit, file_name, fr_am)

    implicit none

    integer, intent(in) :: unit
    character (len=*), intent(in) :: file_name
    type (harm_freq_ampl), dimension (:), intent(in) :: fr_am
    integer (kind=iac) :: num_harmonics, ii
    integer :: io

    num_harmonics = size(fr_am)

    open (unit=unit, file=file_name, action='write', status='replace')

    do ii=0,num_harmonics

      if (ii==0) then

        write(unit=unit, fmt="(9X,A,4X,A,15X,A,16X,A,13X,A)", iostat=io) "# ", " Frequencies ", " Amplitudes ", &
                                                                      " Re[Amplitude] ", " Im[Amplitude] "
        write(unit=unit, fmt=*, iostat=io) ""

      else

        write(unit=unit, fmt="(I10,5X,ES23.16E2,5X,ES23.16E2,5X,ES23.16E2,5X,ES23.16E2)", iostat=io) fr_am(ii)%harmonic_num,  &
                                           fr_am(ii)%frequency, fr_am(ii)%amplitude,fr_am(ii)%real_ampl, fr_am(ii)%imaginary_ampl

      endif

      if (io > 0) then

        print*, "ERROR - subroutine StoreFreqAmp !!!"
        print*, "Some error occured while writing the file ", file_name, " with unit: ", unit, " !!!"

      endif

    enddo

    !#!print*, "The file ", file_name, " with unit: ",unit ," is written normally."

    close (unit=unit)


  end subroutine StoreFreqAmp

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


end module READ_STORE_DATA

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

module NAFF_UV

  use, intrinsic :: iso_c_binding
  use KIND_ACCURACY
  use NUMBERS_CONSTANTS
  use NEEDED_STRUCTURES

  implicit none

contains

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !+
  ! function Naff(input_csignal, sampling_frequency, type_of_inpout_signal, num_freqs, obtain_zero, f_win, a_win) result(freqs_amps)
  !
  ! This function implements the NAFF algorithm for calculating the spectra of periodic data.
  !
  ! See NAFF_UV documentation for details.
  !
  ! Frequencies returned are in units of 2pi.  i.e. freqs ranges from 0 to 1 if the input is complex and between 0 to 0.5 if it is real.
  !
  ! This subroutine will repeat the decomposition loop until all asked elements of freqs and amps are found.
  !
  ! Input:
  !   input_csignal(:) :       -> complex(kind=rac), complex signal for complex data, complex signal with zero the imajinary part for real data
  !   sampling_frequency :     -> 1/sample_spacing, f_n=freq*sampling_frequency and if sampling_frequency<1 we may not be able to retrieve
  !                               the full decimal part of the tune. In contrast if sampling_frequency >=1 there is not such a problem
  !   type_of_inpout_signal :  -> 'c' if it is complex and 'r' if it is real
  !   num_freqs :              -> number of the freqiencies nedded
  !   obtain_zero :            -> if we need to store the dc component
  !   f_win :                  -> the window that scale the data for better determination of the frequency
  !   a_win: 			       -> the window that scale the data for better determination of the amplitude
  !   mod_gs: 				   -> switch for the use of the modified Gram - Schmidt orthonormalization
  !
  ! Output:
  !   freqs_amps(:) :          -> special type that store the harmonic, frequency, amplitude, real and imaginary part of the amplitude
  !                               harmonic=0 for dc / harmonic=1 for the freq. with the largest ampl. and so on
  !-

  function Naff(input_csignal, sampling_frequency, type_of_inpout_signal, &
  				num_freqs, obtain_zero, f_win, a_win, mod_gs) result(freqs_amps)

    implicit none

    complex (kind=rac), dimension(:), intent(in) :: input_csignal
    real (kind=rac), intent(in) :: sampling_frequency
    character (len=*), intent(in) :: type_of_inpout_signal
    integer (kind=iac), intent(in) :: num_freqs, obtain_zero, mod_gs
    character (len=*), intent(in) :: f_win, a_win
    type (harm_freq_ampl), dimension(:), allocatable :: freqs_amps

    complex (kind=rac), dimension(:), allocatable :: dummy_csignal
    complex (kind=rac), dimension(:,:), allocatable :: record
    integer (kind=iac) :: num_of_data, max_num_freqs, ii, jj, time
    real (kind=rac) :: freq, norm_coef
    complex (kind=rac) :: camp


    num_of_data = size(input_csignal)

    max_num_freqs=num_freqs+1

    allocate(dummy_csignal(num_of_data))
    allocate(record(num_of_data,max_num_freqs))
    allocate(freqs_amps(num_freqs+obtain_zero))

    dummy_csignal=input_csignal

    do ii=1,max_num_freqs

      if(ii==1) then
        freq = 0.0_rac
      else
        ! estimate location of spectrum peak using FFTW
        freq = InterpolatedFFTW (dummy_csignal, num_of_data, type_of_inpout_signal, f_win)
        ! refine location of frequency peak
        freq = OptimizeProjection (freq, dummy_csignal, num_of_data, f_win)
      endif
      
      ! construct the harmonic
      record(:,ii) = EigenStateData(freq,num_of_data)

      ! modified Gram - Schmidt orthogonalization
      if (mod_gs==1_iac) then

      	do jj=3,ii
      	  
      	  ! making use of the a_win in modified Gram - Schmidt
          record(:,ii) = record(:,ii)-(ScalarProduct(record(:,ii),ApplyWindow(record(:,jj-1),num_of_data,a_win)) &
                         / sum(conjg(record(:,jj-1))*record(:,jj-1)))*record(:,jj-1)
          
          ! standar modified Gram - Schmidt (without any window)
          !record(:,ii) = record(:,ii)-(ScalarProduct(record(:,ii),record(:,jj-1)) &
          !              / ScalarProduct(record(:,jj-1),record(:,jj-1)))*record(:,jj-1)
      		
      	enddo
      	!normalization
        !record(:,ii) = record(:,ii)/sqrt(conjg(record(:,ii))*record(:,ii))

      endif
      
      
      ! calculate the amplitude with the use of a window and remove the found harmonic
      if (ii>1 .and. type_of_inpout_signal=='r') then
        
        camp = 2.0_rac*ScalarProduct(record(:,ii), ApplyWindow(dummy_csignal, num_of_data, a_win))
        
        ! remove the found harmonic
        dummy_csignal = dummy_csignal - CosData (freq, camp, num_of_data)

      else

        camp = ScalarProduct(record(:,ii), ApplyWindow(dummy_csignal, num_of_data, a_win))

        ! remove the found harmonic
        dummy_csignal = dummy_csignal - camp*record(:,ii)

      endif
      
      
      ! store frequency and amplitude
      if (obtain_zero==1) then
        
        freqs_amps(ii)%harmonic_num=ii-1
        freqs_amps(ii)%frequency=freq*sampling_frequency
        freqs_amps(ii)%amplitude=abs(camp)
        freqs_amps(ii)%real_ampl=real(camp)
        freqs_amps(ii)%imaginary_ampl=aimag(camp)
        
      else

        if (ii==1) then

          cycle

        else
          
          freqs_amps(ii-1)%harmonic_num=ii-1
          freqs_amps(ii-1)%frequency=freq*sampling_frequency
          freqs_amps(ii-1)%amplitude=abs(camp)
          freqs_amps(ii-1)%real_ampl=real(camp)
          freqs_amps(ii-1)%imaginary_ampl=aimag(camp)
          
        endif

      endif
    
    enddo

  end function Naff

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function ScalarProduct (ket_a, ket_b) result (projection)

    implicit none

    complex (kind=rac), dimension(:), intent(in) :: ket_a, ket_b
    complex (kind=rac) :: projection
    complex (kind=rac), dimension(size(ket_a)) :: bra

    bra=conjg(ket_a)
    projection = sum(bra*ket_b)/size(ket_a)

  end function ScalarProduct

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function EigenStateData (fj, data_num) result (wave_data)

    implicit none

    real (kind=rac), intent(in) :: fj
    integer (kind=iac), intent(in) :: data_num
    complex (kind=rac), dimension(data_num) :: wave_data
    integer (kind=iac) :: time

    wave_data = (/ (exp((0.0_rac,1.0_rac)*twopi*fj*time),time=0,data_num-1) /)

  end function EigenStateData

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function CosData (fj, c_amp, data_num) result (cos_wave_data)

    implicit none

    real (kind=rac), intent(in) :: fj
    complex (kind=rac), intent(in) :: c_amp
    integer (kind=iac), intent(in) :: data_num
    real (kind=rac), dimension(data_num) :: cos_wave_data
    real (kind=rac) :: phase
    integer (kind=iac) :: time

    phase = atan2(aimag(c_amp),real(c_amp))

    cos_wave_data = (/ ((abs(c_amp)*cos((twopi*fj*time)+phase)),time=0,data_num-1) /)

  end function CosData

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function ApplyWindow (input_data, num_of_in_data, window_type) result (out_wdata)

    implicit none

    complex (kind=rac), dimension(:), intent(in) :: input_data
    integer (kind=iac), intent(in) :: num_of_in_data
    character (len=*), intent(in) :: window_type
    complex (kind=rac), dimension(num_of_in_data) :: out_wdata

    integer (kind=iac) :: ii, jj, pp
    real (kind=rac) :: cnpp, ss, rr, zz
    real (kind=rac), dimension(num_of_in_data) :: window

    !hann p
    if (window_type(1:3)=='han') then

      read(window_type(5:),'(i10)') pp
      cnpp=1.0_rac
      do jj=1, pp
        cnpp=cnpp*2._rac*jj/(pp+jj)
      enddo

      do ii=1, num_of_in_data

        window(ii) = cnpp*(1.0_rac + cos(twopi*(((ii-1.0_rac)/(num_of_in_data))-0.5_rac)))**pp

      enddo

    !gaussian
    elseif (window_type(1:3)=='gau') then

      read(window_type(5:),*) ss

      do ii=1, num_of_in_data

        rr = (num_of_in_data-1.0_rac)/2.0_rac
        window(ii) = exp(-0.5_rac*((1.0_rac*ii-rr-1.0_rac)/(ss*rr))**2)

      enddo

    ! flat-top -> HFT248D
    elseif (window_type(1:3)=='hft') then

      do ii=1, num_of_in_data

        zz=twopi*((ii-1.0_rac)/(num_of_in_data))

        window(ii) = one - 1.985844164102_rac*cos(zz) + 1.791176438506_rac*cos(two*zz) - 1.282075284005_rac*cos(three*zz) &
                      + 0.667777530266_rac*cos(four*zz) - 0.240160796576_rac*cos(five*zz) + 0.056656381764_rac*cos(six*zz) &
                      - 0.008134974479_rac*cos(seven*zz) + 0.000624544650_rac*cos(eight*zz) - 0.000019808998_rac*cos(nine*zz) &
                      + 0.000000132974_rac*cos(ten*zz)

      enddo

    ! no window
    elseif (window_type(1:3)=='rec') then

      out_wdata = input_data
      return

    else

      print*, "WARNING - ApplyWindow"
      print*, "The fierst three letters of the provided window name, ", window_type," is not supported!"
      print*, "Instead the 'rec' (rectangular window) is applied!"

      out_wdata = input_data
      return

    endif

    out_wdata = input_data * window


  end function ApplyWindow

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ! function maximize_projection
  !
  ! Optimizer that uses an appropriate window (in time domain) and the mnbrak and brent procedures to find a global maximum, which is the frequency that maximizes the projection.


  function OptimizeProjection(f_ini_guess, in_data, number_of_in_data, type_of_fr_window) result (out_freq)

    use BRENT_MNBRAK_METHODS

    implicit none

    real (kind=rac) :: f_ini_guess
    complex (kind=rac), dimension(:), intent(in) :: in_data
    integer (kind=iac), intent(in) :: number_of_in_data
    character (len=*), intent(in) :: type_of_fr_window
    real (kind=rac) :: out_freq

    complex (kind=rac) :: out_camp
    type(brent_mnbrak_class) :: myfunc
    complex (kind=rac), dimension(number_of_in_data) :: in_wdata
    real (kind=rac) :: fr_a, fr_b, fr_c, fun_1, fun_2, fun_3
    real (kind=rac) :: fr_d, fr_e, fr_f, tol, dummy_amp_out_freq


    !apply window for accurate calculation of the freqiencies
    
    if ((f_ini_guess - 0.1_rac/number_of_in_data) <= 0.0_rac) then
      	fr_a = f_ini_guess
        fr_b = f_ini_guess + 0.1_rac/number_of_in_data
    else
    	fr_a = f_ini_guess
        fr_b = f_ini_guess - 0.1_rac/number_of_in_data
    endif

    if (type_of_fr_window(1:3)=='hft') then

      in_wdata = ApplyWindow (in_data, number_of_in_data, 'han_1')

      call myfunc%set_function(MirrorProjection)

      call myfunc%braket_min(fr_a, fr_b, fr_c, fun_1, fun_2, fun_3)

      tol = epsilon (tol)

      fr_d = min(fr_a, fr_b, fr_c)
      fr_e = max(fr_a, fr_b, fr_c)
      fr_f =  fr_a + fr_b + fr_c - fr_d - fr_e
      if (fr_d < 0.0_rac) then
         fr_d = tol
      endif
      if (fr_f < 0.0_rac) then
        fr_f = (fr_e+fr_d)/2.0_rac
      endif

      call myfunc%global_min(fr_d, fr_e, fr_f, tol, out_freq, dummy_amp_out_freq)

      in_wdata = ApplyWindow (in_data, number_of_in_data, type_of_fr_window)

      call myfunc%set_function(Projection)

      fr_d = out_freq - 0.1_rac/number_of_in_data
      fr_e = out_freq + 0.1_rac/number_of_in_data
      fr_f = out_freq

      call myfunc%global_min(fr_d, fr_e, fr_f, tol, out_freq, dummy_amp_out_freq)

    else

      in_wdata = ApplyWindow (in_data, number_of_in_data, type_of_fr_window)

      call myfunc%set_function(MirrorProjection)

      call myfunc%braket_min(fr_a, fr_b, fr_c, fun_1, fun_2, fun_3)

      tol = epsilon (tol)

      fr_d = min(fr_a, fr_b, fr_c)
      fr_e = max(fr_a, fr_b, fr_c)
      fr_f = fr_a + fr_b + fr_c - fr_d - fr_e
      if (fr_d < 0.0_rac) then
         fr_d = tol
      endif
      if (fr_f < 0.0_rac) then
        fr_f = (fr_e+fr_d)/2.0_rac
      endif
      
      call myfunc%global_min(fr_d, fr_e, fr_f, tol, out_freq, dummy_amp_out_freq)

    endif    


    contains

    ! function Projection
    ! Calculates <exp(i*fr*jj)|in_wdata> jj->t
    ! Used only by OptimizeProjection.

    function Projection (me,x) result (f)
      class(brent_mnbrak_class),intent(inout) :: me
      real (kind=rac), intent(in) :: x
      real (kind=rac) :: f
      f = 1.0_rac*abs(ScalarProduct(EigenStateData(x,number_of_in_data),in_wdata))
    end function Projection
    

    function MirrorProjection (me,x) result (f)
      class(brent_mnbrak_class),intent(inout) :: me
      real (kind=rac), intent(in) :: x
      real (kind=rac) :: f
      f = -1.0_rac*abs(ScalarProduct(EigenStateData(x,number_of_in_data),in_wdata))
    end function MirrorProjection

  end function OptimizeProjection

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !  function InterpolatedFFTW
  !
  !  Windows the data and used fftw to find the peak in the spectrum.
  !  The result is interpolated to improve the accuracy.  Hann, Gaussian and Flat-top windowing are
  !  available.


  function InterpolatedFFTW (in_data, number_of_in_data, type_of_in_data, type_of_window) result (out_freq)

    use, intrinsic :: iso_c_binding
    use KIND_ACCURACY

    implicit none

    include 'fftw3.f08'

    complex (kind=rac), dimension(:), intent(in) :: in_data
    integer (kind=iac), intent(in) :: number_of_in_data
    character (len=*), intent(in) :: type_of_in_data
    character (len=*), intent(in) :: type_of_window
    real (kind=rac) :: out_freq

    integer (kind=iac) :: max_ind, cout_size
    complex (kind=rac), dimension(number_of_in_data) :: in_wdata
    real (kind=rac) :: rr, amp_va, amp_va_m, amp_va_p, cm_amp
    real (kind=rac), dimension(:), allocatable :: fftw_amp


    type (C_PTR) :: plan_fw, ptin, ptout
    complex (C_DOUBLE_COMPLEX), pointer :: cin(:), cout(:)
    real(C_DOUBLE), pointer :: rin(:)
    integer (C_INT), dimension(1) :: nn


    if (number_of_in_data<30) then

      Print*, "WARNING - function InterpolatedFFTW !!!"
      print*, "The number of signal data must be more than 30 for an accurate calculation of the frequence !!!"

    endif

    !apply window
    if (type_of_window(1:3)=='hft') then

      in_wdata = ApplyWindow (in_data, number_of_in_data, 'han_1')

    else

      in_wdata = ApplyWindow (in_data, number_of_in_data, type_of_window)

    endif


    nn(1)=int(number_of_in_data,kind=C_INT)

    if (type_of_in_data=='c') then

      cout_size=number_of_in_data

      ptin = fftw_alloc_complex(int(number_of_in_data,C_SIZE_T))
      call c_f_pointer (ptin, cin, [number_of_in_data])
      ptout = fftw_alloc_complex(int(cout_size,C_SIZE_T))
      call c_f_pointer (ptout, cout, [cout_size])


      cin = cmplx(in_wdata, kind=C_DOUBLE_COMPLEX)

      plan_fw = fftw_plan_dft (size(shape(nn)), nn, cin, cout, FFTW_FORWARD, FFTW_ESTIMATE)

      call fftw_execute_dft ( plan_fw, cin, cout)

      call fftw_destroy_plan (plan_fw)

      allocate(fftw_amp(cout_size))


    elseif (type_of_in_data=='r') then

      cout_size=int(number_of_in_data/2_iac + 1_iac, kind=iac)

      ptin = fftw_alloc_real(int(number_of_in_data,C_SIZE_T))
      call c_f_pointer (ptin, rin, [number_of_in_data])
      ptout = fftw_alloc_complex(int(cout_size,C_SIZE_T))
      call c_f_pointer (ptout, cout, [cout_size])


      rin = real(in_wdata, kind=C_DOUBLE)

      plan_fw = fftw_plan_dft_r2c (size(shape(nn)), nn, rin, cout, FFTW_ESTIMATE)

      call fftw_execute_dft_r2c (plan_fw, rin, cout)

      call fftw_destroy_plan (plan_fw)

      allocate(fftw_amp(cout_size))


    else
      print*, "ERROR - FFTW !!!"
      print*, "The type of the input data must be complex <c> or real <r> !!!"
      return
    endif

    fftw_amp(:)=real (abs(cout(:)), kind=rac)

    max_ind = maxloc(fftw_amp(2:), 1) + 1 ! the DC is already removed

    if (fftw_amp(max_ind) == 0.0_rac) then
      print*, "No other freqiencies can be found in the given signal so scanning is terminated !!!"
      return
    endif

    if (max_ind==cout_size) then

      out_freq = 1.0_rac*(max_ind-1_iac)/number_of_in_data

    else

      if (type_of_window(1:3)=='han' .or. type_of_window(1:3)=='hft' .or. type_of_window(1:3)=='rec') then

        ! Parabolic Interpolation (use with hann window)
        amp_va = fftw_amp(max_ind)
        amp_va_m = fftw_amp(max_ind-1_iac)
        amp_va_p = fftw_amp(max_ind+1_iac)
        cm_amp = ((amp_va_p-amp_va_m)/2.0_rac)/(2.0_rac*amp_va - amp_va_p - amp_va_m)
        out_freq = (1.0_rac*(max_ind-1_iac) + cm_amp)/number_of_in_data

      elseif (type_of_window(1:3)=='gau') then

        ! Gaussian Interpolation (use with gaussian window)
        amp_va = log(fftw_amp(max_ind))
        amp_va_m = log(fftw_amp(max_ind-1_iac))
        amp_va_p = log(fftw_amp(max_ind+1_iac))
        cm_amp = ((amp_va_p-amp_va_m)/2.0_rac)/(2.0_rac*amp_va - amp_va_p - amp_va_m)
        out_freq = (1.0_rac*(max_ind-1_iac) + cm_amp)/number_of_in_data

      endif

    endif

    call fftw_free(ptin)
    call fftw_free(ptout)

  end function InterpolatedFFTW

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end module NAFF_UV

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


program naff_uv_main

use KIND_ACCURACY
use NEEDED_STRUCTURES
use READ_STORE_DATA
use NAFF_UV

implicit none

integer :: io,read_skip=1
character (len=:), allocatable :: control_data
character (len=100) :: control_data0, fname_in_data0, win0_f, win0_a, fname_stored_freq_amp0
character (len=1) :: signal_type
character (len=:), allocatable :: fname_in_data, win_f, win_a, fname_stored_freq_amp
character (len=4) :: full_part_signal
integer (kind=iac) :: numb_of_ini_acti_line, numb_of_needed_data, f_p_data_step
real (kind=rac) :: fs
integer (kind=iac) :: number_of_frequensies, dc, mod_gs_orthnorm
type(input_data), dimension(:), allocatable :: signal_for_analyze
complex (kind=rac), dimension(:), allocatable :: input_signal
type(harm_freq_ampl), dimension(:), allocatable :: signal_spectrum
integer (kind=iac) :: numb_lines, numb_active_lines, f_p_numb_active_lines


!read naff_uv_control_file
  open(unit=17,file='naff_uv_control_file.txt',status='old', action='read')

    do

      read(unit=17, fmt=*, iostat=io) control_data0
      control_data = trim(control_data0)

      if (io > 0) then

        print*, "ERROR !!!"
        print*, "Some error occured while reading the naff_uv_control_file.txt !!!"
        rewind(17)
        return

      else if (io == 0 .and. control_data(1:1) /= '#') then

        if (read_skip==1) then
          read (control_data, fmt=*) fname_in_data0
          fname_in_data = trim(fname_in_data0)
        elseif (read_skip==2) then
          read (control_data, fmt=*) signal_type
        elseif (read_skip==3) then
          read (control_data, fmt=*) full_part_signal
        elseif (read_skip==4) then
          read (control_data, fmt=*) numb_of_ini_acti_line
        elseif (read_skip==5) then
          read (control_data, fmt=*) numb_of_needed_data
        elseif (read_skip==6) then
          read (control_data, fmt=*) f_p_data_step
        elseif (read_skip==7) then
          read (control_data, fmt=*) fs
        elseif (read_skip==8) then
          read (control_data, fmt=*) number_of_frequensies
        elseif (read_skip==9) then
          read (control_data, fmt=*) dc
        elseif (read_skip==10) then
          read (control_data, fmt=*) win0_f
          win_f = trim(win0_f)
        elseif (read_skip==11) then
          read (control_data, fmt=*) win0_a
          win_a = trim(win0_a)
        elseif (read_skip==12) then
          read (control_data, fmt=*) mod_gs_orthnorm
        elseif (read_skip==13) then
          read (control_data, fmt=*) fname_stored_freq_amp0
          fname_stored_freq_amp = trim(fname_stored_freq_amp0)
        endif

        read_skip = read_skip + 1

      else if (io == 0 .and. control_data(1:1) == '#') then

        cycle

      else

        !#!print*, "The naff_uv_control_file.txt is read normally."
        exit

      endif

    enddo

  rewind (17)

  close (17)

  !#!print*, '#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#'
  !#!print*, '#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#'
  !#!print*, 'The signal is taken from the file: ', fname_in_data
  !#!print*, '#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#'
  !#!print*, '#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#'

  allocate (signal_spectrum(number_of_frequensies+dc))


  !read data_file
  call GetLines (17, fname_in_data, numb_lines, numb_active_lines)

  if (full_part_signal=='full') then

    numb_of_ini_acti_line = 1
    f_p_numb_active_lines = int((1.0_rac*numb_active_lines/f_p_data_step), kind=iac)
    !#!print*, 'Based on the provided data step: ',f_p_data_step, ' the number of the used data is: ', f_p_numb_active_lines

  elseif (full_part_signal=='part') then

    if ((numb_active_lines-numb_of_ini_acti_line+1) < ((numb_of_needed_data - 1)*f_p_data_step + 1)) then

      f_p_numb_active_lines = int((1.0_rac*(numb_active_lines-numb_of_ini_acti_line+1))/f_p_data_step, kind=iac)
      print*, "WARNING !!!"
      print*, "naff_main: The number of the asked data are more than the avelable ones !!!"
      print*, "Thus, the maximum number of data with initial row: ", numb_of_ini_acti_line, " and data step: ", &
              f_p_data_step, " is used and it is equal to: ", f_p_numb_active_lines

    else

      f_p_numb_active_lines = numb_of_needed_data

    endif

  endif

  allocate (signal_for_analyze(f_p_numb_active_lines))
  allocate (input_signal(f_p_numb_active_lines))

  signal_for_analyze = GetData (17, fname_in_data, numb_of_ini_acti_line, f_p_numb_active_lines, f_p_data_step, signal_type)


  !start naff
  input_signal = cmplx(signal_for_analyze%real_part, signal_for_analyze%imaginary_part, kind=rac)

  signal_spectrum = Naff(input_signal, fs, signal_type, number_of_frequensies, dc, win_f, win_a, mod_gs_orthnorm)


  !store asked freqiencies and amplitude
  call StoreFreqAmp (17, fname_stored_freq_amp, signal_spectrum)

end program naff_uv_main
