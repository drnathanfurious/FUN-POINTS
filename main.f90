! FUN-POINTS
! 
! Solution to this Erdos thing.  You could win $50!

program main

  ! this just makes things pretty and keeps all the derived types in one place
  use points_module
  use optimization_module

 implicit none


 integer,parameter :: number_of_points = 6
 integer,parameter :: dimensions = 2


 !!!!!
 ! a flattened version of all the positions of the points
 ! fortran matrices are just arrays... but in >>> COLUMN-MAJOR ORDER <<<
 real, target :: points(number_of_points,dimensions)  ! let's just stick with arrays...
 ! using the idea of an upper triangular adjaceny matrix to track
 ! distances between points...
 real,target :: adjacency_matrix (number_of_points,number_of_points)

 type(opt_function_type) :: opt_function

 opt_function%ires=0

 ! randomly select points in the plane
 points = InitializePoints(number_of_points,dimensions)

 ! initialize the optimization routine
 call InitOptimizationFunction (opt_function)

 ! assign the points and adjacency matrix to the opt_function structure
 opt_function%f_data%points => points
 opt_function%f_data%adjacency_matrix => adjacency_matrix



 !!!!!! start optimization
 call nlo_optimize(opt_function%ires, opt_function%opt, points, opt_function%minf)

 if (opt_function%ires.lt.0) then 
   write(*,*) "nlopt failed! This reflects poorly on you."
 else
   write(*,*) "found min at "
   write (*,*) points
   write(*,*) "minimum ", opt_function%minf
 end if

 !!!!!! clean up
 call nlo_destroy(opt_function%opt)
 adjacency_matrix = UpdateAdjacencyMatrix(points)
 !call PrintResults(points)
 !i = PrintDistances(adjacency_matrix)


end program main

