program finite_diff
implicit none      

    ! type declaration statements
    integer::t_choice,x_choice,Nx,Nt,Final_time,error 
    real :: X_min,X_max,h,dt,dt_over_h
    logical :: impl
    real, allocatable :: X(:),u1H0(:),u2H0(:),u3H0(:),PH0(:),next_u1(:),next_u2(:),next_u3(:)

    ! time discretization: 
    ! du/dt(t^n) =  
    ! 1:  ( u^n - u^n-1 )/dt        (impl) 
    ! 2:  ( u^n+1 - u^n )/dt        (expl) 
    t_choice = 2 
    
    ! space discretization: 
    !du/dx(x_i) =  
    !1:  ( u_i - u_i-1 )/dx 
    !2:  ( u_i+1 - u_i )/dx 
    !3:  ( u_i+1 - u_i-1 )/2dx
    x_choice = 3
    
    ! spatial domain and mesh
    X_min = 0
    X_max = 1
    Nx = 100
    h = 1./Nx
    allocate(X(Nx),stat=error)
    if (error.ne.0) then
        print*,'error: could not allocate memory for array, Nx=',Nx
        stop
    endif

    ! temporal domain
    Final_time = 1
    Nt = 120
    dt = 1./Nt ! Final_time * 1./Nt
    dt_over_h = 1.0*dt/h
    
    allocate(u1H0(Nx),stat=error)
    if (error.ne.0) then
        print*,'error: could not allocate memory for array, Nx=',Nx
        stop
    endif
    
    allocate(PH0(Nx),stat=error)
    if (error.ne.0) then
        print*,'error: could not allocate memory for array, Nx=',Nx
        stop
    endif


    ! executable statements  
    print*,"Some results Nt",Nt
    
    ! deallocation and end of program
    deallocate(X,stat=error)
    if (error.ne.0) then
        print*,'error in deallocating array'
    endif
    
    deallocate(u1H0,stat=error)
    if (error.ne.0) then
        print*,'error in deallocating array'
    endif
    
    deallocate(PH0,stat=error)
    if (error.ne.0) then
        print*,'error in deallocating array'
    endif
end program finite_diff
