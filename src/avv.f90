      module avv_mod
  
      contains

      function avv(ax1,ay1,az1,ax2,ay2,az2,ax3,ay3,az3,ax4,ay4,az4)
	use params
	real, value :: ax1,ay1,az1,ax2,ay2,az2,ax3,ay3,az3,ax4,ay4,az4
        real :: ax12,ay12,az12,ax34,ay34,az34,avv
!$acc routine seq
!$OMP declare target 
      ax12=ax2-ax1
      ay12=ay2-ay1
      az12=az2-az1
      ax34=ax4-ax3
      ay34=ay4-ay3
      az34=az4-az3
      avv=(ax12*ax34+ay12*ay34+az12*az34)/(0.000001+
     $     sqrt((ax12**2+ay12**2+az12**2)*(ax34**2+ay34**2+az34**2)))

c^^^^^^^^^^^^^ v.v finished ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end function avv
      end module avv_mod
