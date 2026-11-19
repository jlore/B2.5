      module b2mod_output_ornl
      implicit none
c ORNL verbosity for optional text diagnostics: 0 is minimum, 10 maximum.
c The driver reads this value before entering the parallel solver loops.
      integer, save :: output_level_ornl = 0
      end module b2mod_output_ornl
