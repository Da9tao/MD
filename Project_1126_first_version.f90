    !  MD_practice_02.f90
    !
    !  FUNCTIONS:
    !  MD_practice_02 - Entry point of console application.
    !

    !****************************************************************************
    !
    !  PROGRAM: MD_practice_02
    !
    !  PURPOSE:  Entry point for the console application.
    !
    !****************************************************************************

    module MD
    implicit none
    real(kind=8), parameter :: N_av=6.02214076d23, D_0=0.1848d0*4.184d3/N_av, R_0=3.5532d-10, m_0=16d-3/N_av
    real(kind=8), parameter :: kb=1.380649d-23, T=298d0
    real(kind=8), parameter :: E_0=1.5d0*kb*T, D_null=D_0/E_0, p_0=E_0/(R_0*R_0*R_0)
    real(kind=8), parameter :: t_0=SQRT(m_0*R_0*R_0/E_0)
    real(kind=8), parameter :: pi=ACOS(-1d0)
    real(kind=8), parameter :: e_charge = 1.602176487d-19, epsilon_vacuum = 8.854187817d-12, electric_energy_null = e_charge * e_charge / (4d0 * pi * epsilon_vacuum * R_0 * E_0)
    real(kind=8), parameter :: bond_length_eq_null = 1d-10/R_0, k_bond = 500d0*4.184d3/N_av*1d20, k_bond_null=k_bond*R_0*R_0/E_0
    real(kind=8), parameter :: cos_theta_eq = -1d0/3d0, sin_theta_eq_square = 8d0/9d0, k_angle = 120d0*4.184d3/N_av, k_angle_null = k_angle/E_0
    real(kind=8), parameter :: Q_nh=40d0*4.184d3/N_av*1d-24, Q_nh_null=Q_nh/(E_0*t_0*t_0)
    real(kind=8), parameter :: E_unit=N_av/4.184d3
    integer, parameter :: seed(2)=(/123, 789/)
    type atom
        real(8) :: D0
        real(8) :: R0
        real(8) :: mass
        real(8) :: charge
    end type

    contains

    subroutine init_position_square(length, Np, Mat)
    implicit none
    real(8), intent(in) :: length
    integer, intent(in) :: Np
    real(8), dimension(3, Np), intent(out) :: Mat
    integer :: i, j, k, counter
    real(8) :: spacing
    integer :: number_edge

    number_edge=INT(real(Np, kind=8)**(1d0/3d0))+1
    spacing=length/real(number_edge, kind=8)
    counter=1
    do i=1, number_edge
        do j=1, number_edge
            do k=1, number_edge
                Mat(1,counter)=real(i, kind=8)*spacing
                Mat(2,counter)=real(j, kind=8)*spacing
                Mat(3,counter)=real(k, kind=8)*spacing
                counter=counter+1
                if (counter>Np) then
                    Mat=Mat-0.5d0*spacing
                    return
                end if
            end do
        end do
    end do
    end subroutine

    subroutine init_position_square_H2O(length, Np, Mat, distance)
    implicit none
    real(8), intent(in) :: length, distance
    integer, intent(in) :: Np
    real(8), dimension(3, Np), intent(out) :: Mat
    integer :: i, j, k, counter
    real(8) :: spacing
    integer :: number_edge
    real(8) :: x_comp, y_comp

    number_edge = INT(NINT(Np/3d0)**(1d0/3d0))+1
    spacing = length / real(number_edge, kind=8)
    counter = 1
    x_comp = SQRT(2d0/3d0)
    y_comp = SQRT(1d0/3d0)

    do k=1, number_edge
        do j=1, number_edge
            do i=1, number_edge
                Mat(1,counter:(counter+2)) = real(i, kind=8)*spacing
                Mat(2,counter:(counter+2)) = real(j, kind=8)*spacing
                Mat(3,counter:(counter+2)) = real(k, kind=8)*spacing

                Mat(1, (counter+1)) = Mat(1, (counter+1)) - x_comp * distance
                Mat(2, (counter+1)) = Mat(2, (counter+1)) + y_comp * distance

                Mat(1, (counter+2)) = Mat(1, (counter+2)) + x_comp * distance
                Mat(2, (counter+2)) = Mat(2, (counter+2)) + y_comp * distance

                counter=counter+3
                if (counter>Np) then
                    Mat = Mat - 0.5d0 * spacing
                    return
                end if
            end do
        end do
    end do

    end subroutine init_position_square_H2O


    subroutine init_velocity(Np, Mat)
    implicit none
    integer, intent(in) :: Np
    real(8), dimension(4, Np), intent(out) :: Mat
    integer counter, i, seed(2)
    real(8) :: temp(3)=0d0, temp_sum=0d0, temp_kinetic=0d0, center_velocity(3)=0d0, center_kinetic

    call RANDOM_NUMBER(Mat(1:3,:))
    Mat(1:3,:)=2d0*Mat(1:3,:)-1d0
    temp_sum=sum(0.5d0*Mat(1:3,:)*Mat(1:3,:))
    !print *, temp_sum
    center_velocity=SUM(Mat(1:3,:),2)/real(Np, kind=8)
    center_kinetic=0.5d0*Np*SUM(center_velocity*center_velocity)
    temp_sum=temp_sum-center_kinetic
    do counter=1, 3
        Mat(counter,:)=Mat(counter,:)-center_velocity(counter)
    end do
    Mat(4,:)=SUM(0.5d0*Mat(1:3,:)*Mat(1:3,:),1)
    !print *, SUM( Mat(1:3,:),2)
    Mat(1:3,:)=Mat(1:3,:)*sqrt(real(Np, kind=8)/temp_sum)
    Mat(4,:)=Mat(4,:)*(real(Np, kind=8)/temp_sum)
    !print *, SUM(Mat(4,:))/real(Np, kind=8)
    !print *, SUM(Mat(1:3,:),2)
    end subroutine

    subroutine init_velocity_H2O(Np, mass, Mat)
    implicit none
    integer, intent(in) :: Np
    real(8), intent(in) :: mass(Np)
    real(8), dimension(4, Np), intent(out) :: Mat
    real(8), allocatable :: Mat_temp(:,:)
    integer :: counter, i, N_mol
    real(8) :: temp_sum, temp_kinetic=0d0, center_velocity(3)=0d0, center_kinetic, mass_total

    temp_sum = 0d0
    N_mol = Np / 3
    allocate(Mat_temp(3, N_mol))
    call RANDOM_NUMBER(Mat_temp(1:3,:))
    Mat_temp(1:3,:) = 2d0 * Mat_temp(1:3,:) - 1d0
    do counter = 1, N_mol
        do i = 1, 3
        Mat(1:3, (3 * (counter - 1) + i)) = Mat_temp(1:3, counter) 
        end do       
    end do
    
    mass_total = SUM(mass)
    !call RANDOM_NUMBER(Mat(1:3,:))
    !Mat(1:3,:) = 2d0 * Mat(1:3,:) - 1d0
    do counter = 1, 3
        !Mat(counter,:) = Mat(counter,:) * 1d0 / SQRT(mass)
        temp_sum = temp_sum + 0.5d0 * SUM(mass * Mat(counter,:) * Mat(counter,:))
        center_velocity(counter) = SUM(mass * Mat(counter,:)) / mass_total
    end do
    !print *, temp_sum
    center_kinetic = 0.5d0 * mass_total * SUM(center_velocity * center_velocity)
    temp_sum = temp_sum - center_kinetic
    do counter=1, 3
        Mat(counter,:) = Mat(counter,:) - center_velocity(counter)
    end do
    Mat(4,:) = 0.5d0 * mass * SUM(Mat(1:3,:) * Mat(1:3,:), 1)
    !print *, SUM( Mat(1:3,:),2)
    Mat(1:3,:) = Mat(1:3,:) * sqrt(real(Np, kind=8) / temp_sum)
    Mat(4,:) = Mat(4,:) * (real(Np, kind=8) / temp_sum)
    !print *, Mat(4,:)
    !print *, SUM(Mat(4,:))/real(Np, kind=8)
    !Mat(1,:) = Mat(1,:) * mass
    !Mat(2,:) = Mat(2,:) * mass
    !Mat(3,:) = Mat(3,:) * mass
    !print *, SUM(Mat(1:3,:),2)
    deallocate(Mat_temp)
    end subroutine init_velocity_H2O

    subroutine calc_potential(D, Np, x, U, force_m)
    implicit none
    real(8), intent(in) :: D
    integer, intent(in) :: Np
    real(8), intent(in) :: x(3,Np)
    real(8), intent(out) :: U, force_m(3, Np)
    integer :: i, j
    real(8) :: temp=0d0, temp_reciprocal=0d0, dx, dy, dz, dr(3), force_temp, temp_sixth

    U=0d0
    force_m=0d0
    dx=0d0
    dy=0d0
    dz=0d0

    do i=1, Np-1
        do j=i+1, Np
            dx=x(1,i)-x(1,j)
            dy=x(2,i)-x(2,j)
            dz=x(3,i)-x(3,j)
            temp=dx*dx+dy*dy+dz*dz
            temp_reciprocal=(1d0/temp)
            temp_sixth=temp_reciprocal*temp_reciprocal*temp_reciprocal
            dr=(/dx, dy, dz/)

            U=U+D*(temp_sixth*temp_sixth-2*temp_sixth)
            force_temp=12d0*D*(temp_sixth*temp_sixth-temp_sixth)*temp_reciprocal
            dr=force_temp*dr
            force_m(:,i)=force_m(:,i)+dr
            force_m(:,j)=force_m(:,j)-dr

        end do
    end do

    end subroutine

    subroutine calc_potential_period(D, Np, cell_length, r_cut, x, U, force_m, pressure)
    implicit none
    real(8), intent(in) :: D, cell_length, r_cut
    integer, intent(in) :: Np
    real(8), intent(in) :: x(3,Np)
    real(8), intent(out) :: U, force_m(3, Np), pressure
    integer :: i, j
    real(8) :: temp=0d0, temp_reciprocal=0d0, temp_sixth=0d0, dx(3), force_temp(3), temp_force

    U=0d0
    force_m=0d0
    dx=0d0
    pressure=0d0

    do i=1, Np-1
        do j=i+1, Np
            dx=x(1:3,i)-x(1:3,j)
            dx=dx-real(NINT(dx/cell_length), kind=8)*cell_length
            temp=sum(dx*dx)
            if (temp > (r_cut*r_cut)) cycle
            temp_reciprocal=(1d0/temp)
            temp_sixth=temp_reciprocal*temp_reciprocal*temp_reciprocal

            U = U + D * (temp_sixth*temp_sixth - 2d0*temp_sixth)
            temp_force = 12d0 * D * temp_reciprocal * (temp_sixth*temp_sixth-temp_sixth)
            !pressure=pressure+SUM(force_temp*dx*dx)
            force_temp=temp_force*dx
            force_m(:,i)=force_m(:,i)+force_temp
            force_m(:,j)=force_m(:,j)-force_temp
            pressure = pressure + SUM(force_temp * dx)

        end do
    end do

    end subroutine calc_potential_period

    subroutine calc_potential_period_LJ(Np, cell_length, r_cut, type_matrix, atom_type, x, U, force_m, pressure, pair_remaining, number_other)
    implicit none
    real(8), intent(in) :: cell_length, r_cut
    integer, intent(in) :: Np, type_matrix(*)
    real(8), intent(in) :: x(3, Np)
    real(8), intent(out) :: U, force_m(3, Np), pressure
    integer, intent(in) :: number_other, pair_remaining(2, number_other)
    integer :: i, j, atom_index(2), type_index(2)
    real(8) :: temp=0d0, temp_reciprocal=0d0, temp_sixth=0d0, dx(3), force_temp(3), D0_temp, R0_temp
    type(atom), intent(in) :: atom_type(*)

    U=0d0
    force_m=0d0
    !dx=0d0
    pressure=0d0

    do i=1, number_other
        atom_index = pair_remaining(:, i)
        type_index = type_matrix(atom_index)
        if (type_index(1) == type_index(2)) then
            D0_temp = atom_type(type_index(1))%D0
            R0_temp = atom_type(type_index(1))%R0
        else
            D0_temp = SQRT(atom_type(type_index(1))%D0 * atom_type(type_index(2))%D0)
            R0_temp = (atom_type(type_index(1))%R0 + atom_type(type_index(2))%R0) / 2d0
        end if
        dx = x(:, atom_index(1)) - x(:, atom_index(2))
        dx = dx - real(NINT(dx / cell_length), kind=8) * cell_length
        temp = SUM(dx*dx)
        if (temp > (r_cut*r_cut*R0_temp*R0_temp)) cycle
        temp_reciprocal = R0_temp * R0_temp / temp
        temp_sixth = temp_reciprocal*temp_reciprocal*temp_reciprocal

        U = U + D0_temp * (temp_sixth*temp_sixth - 2d0*temp_sixth)
        !temp_force = 12d0 * D0_temp / temp * (temp_sixth*temp_sixth - temp_sixth)
        !pressure=pressure+SUM(force_temp*dx*dx)
        force_temp = 12d0 * D0_temp / temp * (temp_sixth*temp_sixth - temp_sixth) * dx
        force_m(:, atom_index(1)) = force_m(:, atom_index(1)) + force_temp
        force_m(:, atom_index(2)) = force_m(:, atom_index(2)) - force_temp
        pressure = pressure + SUM(force_temp * dx)
    end do
    end subroutine calc_potential_period_LJ


    subroutine calc_Coulomb(Np, cell_length, r_cut, charge, x, U, force_m, pressure, pair_remaining, number_other, alpha, k_matrix, k_counts, pair_bonding, number_bonds, pair_angle, number_angle)
    implicit none
    real(8), intent(in) :: cell_length, r_cut, alpha
    integer, intent(in) :: Np, k_counts, number_bonds, number_angle       ! Np denotes the number of atoms in the system.
    real(8), intent(in) :: x(3,Np), charge(Np), k_matrix(3, k_counts)
    real(8), intent(out) :: U, force_m(3, Np), pressure
    integer, intent(in) :: number_other, pair_remaining(2, number_other), pair_bonding(2, number_bonds), pair_angle(2, number_angle)
    integer :: i, j, atom_index(2), type_index(2)
    real(8) :: U_first_part=0d0, dx(3), force_temp(3, Np), k_temp(3), k_temp_square, rho_k_square, pressure_temp, charge_product, pressure_temp_first
    complex(8) :: rho_k, temp_force(Np), temp_force_extended(3, Np)
    real(8) :: coefficient, distance, distance_square, alpha_root, U_third_part, force_third_part(3, Np), temp_force_third_part(3), alpha_pi_root, temp_factor
    real(8) :: center(3)
    
    force_temp = 0d0
    alpha_root = SQRT(alpha)
    alpha_pi_root =  SQRT(alpha / pi)
    U_first_part = 0d0
    pressure_temp_first = 0d0
    center = 0.5d0 * cell_length

    do i = 1, k_counts
        k_temp = k_matrix(:, i)
        k_temp_square = DOT_PRODUCT(k_temp, k_temp)
        !rho_k = SUM(charge * EXP((0d0, 1d0) * (k_temp(1) * x(1,:) + k_temp(2) * x(2,:) + k_temp(3) * x(3,:))))
        temp_force = charge * EXP((0d0, 1d0) * MATMUL(k_temp, x))  ! Complex
        rho_k = SUM(temp_force)       ! Faster than the previous line.
        rho_k_square = real(rho_k * DCONJG(rho_k), 8)
        U_first_part = U_first_part + 1d0 / k_temp_square * rho_k_square * EXP(- k_temp_square / (4d0 * alpha))

        !temp_force = charge * EXP((0d0, 1d0) * MATMUL(k_temp, x))  ! Complex

        do j = 1, 3
            temp_force_extended(j, :) = temp_force * DCONJG(rho_k) * (0d0, 1d0) * k_temp(j)
        end do
        temp_force_extended = temp_force_extended + DCONJG(temp_force_extended)

        force_temp = force_temp + real(temp_force_extended, 8) / k_temp_square * EXP(- k_temp_square / (4d0 * alpha))
    end do
    
    coefficient = 2d0 * pi * electric_energy_null / (cell_length * cell_length * cell_length)
    U_first_part = U_first_part * coefficient
    force_temp = force_temp * (-coefficient)
    do i = 1, Np
        dx = x(:, i) - center
        pressure_temp_first = pressure_temp_first + DOT_PRODUCT(dx, force_temp(:, i))
    end do
    
    !print *, "U_coumlomb = ", U_first_part
    !print "(3(G0,1X))", force_temp

    ! Second part
    U_first_part = U_first_part - alpha_pi_root * electric_energy_null * SUM(charge * charge)
    !print *, "U_coumlomb = ", U_first_part

    ! Third part
    U_third_part = 0d0
    force_third_part = 0d0
    pressure_temp = 0d0
    do i=1, number_other
        atom_index = pair_remaining(:, i)
        dx = x(1:3, atom_index(1)) - x(1:3, atom_index(2))
        dx = dx - real(NINT(dx/cell_length), kind=8) * cell_length
        distance_square = SUM(dx*dx)
        if (distance_square > (r_cut*r_cut)) cycle
        distance = SQRT(distance_square)
        charge_product = charge(atom_index(1)) * charge(atom_index(2))
        U_third_part = U_third_part + charge_product * ERFC(alpha_root * distance) / distance
        !print *, distance
        temp_factor = 2d0 * alpha_pi_root * EXP(-alpha * distance_square) + ERFC(alpha_root * distance) / distance
        temp_force_third_part(:) = charge_product / distance_square * temp_factor * dx
        force_third_part(:, atom_index(1)) = force_third_part(:, atom_index(1)) + temp_force_third_part
        force_third_part(:, atom_index(2)) = force_third_part(:, atom_index(2)) - temp_force_third_part
        pressure_temp = pressure_temp + charge_product * temp_factor
    end do

    do i=1, number_bonds
        atom_index = pair_bonding(:, i)
        dx = x(1:3, atom_index(1)) - x(1:3, atom_index(2))
        dx = dx - real(NINT(dx/cell_length), kind=8) * cell_length
        distance_square = SUM(dx*dx)
        if (distance_square > (r_cut*r_cut)) cycle
        distance = SQRT(distance_square)
        charge_product = charge(atom_index(1)) * charge(atom_index(2))
        U_third_part = U_third_part - charge_product * ERF(alpha_root * distance) / distance
        !print *, distance
        temp_factor = 2d0 * alpha_pi_root * EXP(-alpha * distance_square) - ERF(alpha_root * distance) / distance
        temp_force_third_part(:) = charge_product / distance_square * temp_factor * dx
        force_third_part(:, atom_index(1)) = force_third_part(:, atom_index(1)) + temp_force_third_part
        force_third_part(:, atom_index(2)) = force_third_part(:, atom_index(2)) - temp_force_third_part
        pressure_temp = pressure_temp + charge_product * temp_factor
    end do

    do i=1, number_angle
        atom_index = pair_angle(:, i)
        dx = x(1:3, atom_index(1)) - x(1:3, atom_index(2))
        dx = dx - real(NINT(dx/cell_length), kind=8) * cell_length
        distance_square = SUM(dx*dx)
        if (distance_square > (r_cut*r_cut)) cycle
        distance = SQRT(distance_square)
        charge_product = charge(atom_index(1)) * charge(atom_index(2))
        U_third_part = U_third_part - charge_product * ERF(alpha_root * distance) / distance
        !print *, distance
        temp_factor = 2d0 * alpha_pi_root * EXP(-alpha * distance_square) - ERF(alpha_root * distance) / distance
        temp_force_third_part(:) = charge_product / distance_square * temp_factor * dx
        force_third_part(:, atom_index(1)) = force_third_part(:, atom_index(1)) + temp_force_third_part
        force_third_part(:, atom_index(2)) = force_third_part(:, atom_index(2)) - temp_force_third_part
        pressure_temp = pressure_temp + charge_product * temp_factor
    end do
    U_third_part = U_third_part * electric_energy_null
    force_third_part = force_third_part * electric_energy_null
    pressure_temp = pressure_temp * electric_energy_null

    U = U + U_first_part + U_third_part
    force_m = force_m + force_temp +  force_third_part
    pressure = pressure + pressure_temp_first + pressure_temp
    !print *, "Third part"
    !print *, U_third_part
    !print *, force_third_part
    !print *, "U_third = ",  U_third_part
    end subroutine calc_Coulomb

    subroutine calc_bond(Np, cell_length, x, U, force_m, pressure, pair_bonding, number_bonds)
    implicit none
    real(8), intent(in) :: cell_length
    integer, intent(in) :: Np
    real(8), intent(in) :: x(3,Np)
    real(8), intent(out) :: U, force_m(3, Np), pressure
    integer, intent(in) :: number_bonds, pair_bonding(2, number_bonds)
    integer :: i, atom_index(2)
    real(8) :: distance=0d0, dx(3), force_temp(3)

    !U=0d0
    !force_m=0d0
    !dx=0d0
    !pressure=0d0

    do i=1, number_bonds
        atom_index = pair_bonding(:, i)
        dx = x(:, atom_index(1)) - x(:, atom_index(2))
        dx = dx - real(NINT(dx/cell_length), kind=8) * cell_length
        distance = sqrt(sum(dx*dx))

        U = U + 1d0/2d0 * k_bond_null * (distance - bond_length_eq_null) * (distance - bond_length_eq_null)
        force_temp = k_bond_null * (bond_length_eq_null / distance - 1d0) * dx
        force_m(:, atom_index(1)) = force_m(:, atom_index(1)) + force_temp
        force_m(:, atom_index(2)) = force_m(:, atom_index(2)) - force_temp
        pressure = pressure + SUM(force_temp * dx)
    end do
    end subroutine calc_bond

    subroutine calc_angle(Np, cell_length, x, U, force_m, pressure, pair_angle, number_angle)
    implicit none
    real(8), intent(in) :: cell_length
    integer, intent(in) :: Np
    real(8), intent(in) :: x(3,Np)
    real(8), intent(out) :: U, force_m(3, Np), pressure
    integer, intent(in) :: number_angle, pair_angle(3, number_angle)
    integer :: i, j, atom_index(3)
    real(8) :: dx(3, 3)=0d0, distance(2), force_temp(3,3), cos_theta, distance_product, pressure_temp(3)
    !U=0d0
    !force_m=0d0
    !dx=0d0
    !pressure=0d0
    do i=1, number_angle
        atom_index = pair_angle(:, i)
        dx(1:3, 1) = x(1:3, atom_index(1)) - x(1:3, atom_index(3))
        dx(1:3, 2) = x(1:3, atom_index(2)) - x(1:3, atom_index(3))
        dx(1:3, 3) = x(1:3, atom_index(1)) - x(1:3, atom_index(2))
        dx = dx - real(NINT(dx/cell_length), kind=8) * cell_length
        distance = SQRT(SUM(dx(:, 1:2) * dx(:, 1:2), 1))
        distance_product = PRODUCT(distance)
        cos_theta = SUM(dx(:, 1) * dx(:, 2)) / distance_product

        U = U + 1d0/2d0 * k_angle_null / sin_theta_eq_square * (cos_theta - cos_theta_eq) * (cos_theta - cos_theta_eq)
        force_temp(:, 1) = - k_angle_null / sin_theta_eq_square * (cos_theta - cos_theta_eq) * (1d0 / distance_product * dx(:,2) - cos_theta / (distance(1)*distance(1)) * dx(:,1))
        force_temp(:, 2) = - k_angle_null / sin_theta_eq_square * (cos_theta - cos_theta_eq) * (1d0 / distance_product * dx(:,1) - cos_theta / (distance(2)*distance(2)) * dx(:,2))
        force_temp(:, 3) = - force_temp(:, 1) - force_temp(:, 2)
        do j = 1, 3
            force_m(:, atom_index(j)) = force_m(:, atom_index(j)) + force_temp(:, j)
        end do
        !pressure = pressure + SUM(force_temp * dx(:, 3))
        pressure_temp(1) = (cos_theta / (distance(1)*distance(1)) - 1d0 / distance_product) * SUM(dx(:, 1) * dx(:, 1))
        pressure_temp(2) = (cos_theta / (distance(2)*distance(2)) - 1d0 / distance_product) * SUM(dx(:, 2) * dx(:, 2))
        pressure_temp(3) = (1d0 / distance_product) * SUM(dx(:, 3) * dx(:, 3))
        pressure = pressure + k_angle_null / sin_theta_eq_square * (cos_theta - cos_theta_eq) * SUM(pressure_temp)
        !print *, force_temp(:, 1)
        !print *, force_temp(:, 2)

    end do
    end subroutine calc_angle

    subroutine pairing(number_H2O, pair_bonding, pair_angle, pair_remaining, molecules_comp)
    integer, intent(in) :: number_H2O
    integer:: number_atoms, number_bonds, number_angle, number_other, number_interaction
    integer, allocatable, intent(out) :: pair_bonding(:,:), pair_angle(:,:), pair_remaining(:,:), molecules_comp(:,:)
    integer :: i, j, k, counter, temp1, temp2, temp_pair(2)
    LOGICAL :: if_unused

    number_atoms = number_H2O * 3
    number_bonds = number_H2O * 2
    number_angle = number_H2O
    number_interaction = number_atoms * (number_atoms-1) / 2
    number_other = number_interaction - number_bonds - number_angle

    allocate(pair_bonding(2, number_bonds))
    allocate(pair_angle(3, number_angle))
    allocate(pair_remaining(2, number_other))
    allocate(molecules_comp(3, number_H2O))
    counter = 1
    pair_bonding = 0
    pair_angle = 0
    pair_remaining = 0
    do i = 1, number_H2O
        pair_bonding(:, counter) = (/3*(i-1) + 1, 3*(i-1) + 2/)
        counter = counter + 1
        pair_bonding(:, counter) = (/3*(i-1) + 1, 3*(i-1) + 3/)
        counter = counter + 1
    end do

    counter = 1

    do i = 1, number_H2O
        pair_angle(:, counter) = (/3*(i-1) + 2, 3*(i-1) + 3, 3*(i-1) + 1/)
        counter = counter + 1
    end do

    counter = 1
    do i = 1, number_atoms-1
        do j = i+1, number_atoms
            if_unused = .true.
            temp_pair=(/i, j/)
            do k = 1, SIZE(pair_bonding, 2)
                if (all(temp_pair == pair_bonding(:, k))) then
                    if_unused = .false.
                    exit
                end if
            end do
            do k = 1, SIZE(pair_angle, 2)
                if (all(temp_pair == pair_angle(1:2, k))) then
                    if_unused = .false.
                    exit
                end if
            end do
            if (if_unused) then
                pair_remaining(:, counter) = temp_pair
                counter = counter + 1
            end if
        end do
    end do

    do i = 1, number_H2O
        molecules_comp(:, i) = (/3*(i-1)+1, 3*(i-1)+2, 3*(i-1)+3/)
    end do

    end subroutine pairing

    subroutine find_points_3d(radius, output, k_counts)
    implicit none
    real(8), intent(in) :: radius
    real(8), allocatable, intent(out) :: output(:,:)
    integer, intent(out) :: k_counts
    real(8), allocatable :: temp_output(:,:)
    real(8), parameter :: pi=ACOS(-1d0)
    integer :: number_estimated, x_max, i, j, k, y_max, z_max, counter

    x_max = INT(radius)
    number_estimated = INT(4d0/3d0 * pi * radius**(3d0) * 1.2d0) + 3
    print *, "Esimated number = ", number_estimated
    allocate(output(3, number_estimated))
    output = 0d0
    counter = 0

    do i = x_max, -x_max, -1
        y_max = int(SQRT(radius**2d0 - real(i, 8)**(2d0)))
        do j = y_max, -y_max, -1
            z_max = int(SQRT(radius**2d0 - real(i, 8)**(2d0) -real(j, 8)**(2d0)))
            do k = z_max, -z_max, -1
                if ((i == 0) .AND. (j == 0) .AND. (k == 0)) cycle
                counter = counter + 1
                output(:, counter) = (/real(i, 8), real(j, 8), real(k, 8)/)
            end do
        end do
    end do

    print *, "Lattice number = ", counter
    allocate(temp_output(3,counter))
    temp_output = output(:, 1:counter)
    deallocate(output)
    allocate(output(3,counter))
    output = temp_output
    k_counts = counter
    deallocate(temp_output)

    end subroutine find_points_3d

    end module MD



    program MD_practice_02
    use MD
    implicit none

    integer :: N_molecules=60, N_atoms                              ! Only water is considered in this case.
    real(8) :: cell_length=7.1d0 , r_cut=3d0, dt=1d-3               ! r_cut means the coefficient.
    integer :: N_steps
    real(8), allocatable :: position(:,:)
    real(8), allocatable :: velocity(:,:)
    real(8), allocatable :: force(:,:), mass(:), mass_matrix(:,:), k_matrix(:,:), charge(:)
    real(8) :: U_potential, U_kinetic, U_total, pressure, pressure_k, rho_null(2), U_correction(2), r_cut_third, pressure_correction, T_mean, volume
    real(8) :: R0_temp, D0_temp, mass_total
    real(8) :: alpha, k_radius
    real(8), dimension(3) :: xyz, center_velocity
    integer :: number_bonds, number_angle, number_other, number_interaction, k_counts
    integer, allocatable :: pair_bonding(:,:), pair_angle(:,:), pair_remaining(:,:), molecules_comp(:,:), type_matrix(:)
    type(atom) :: atom_type(2)

    integer :: counter, i, j, k
    character(len=3) :: residue_name="SOL"
    character(len=2) :: atom_name(3) = (/"O1", "H2", "H3"/)
    real(8) :: start_time, end_time, elapsed_time
    real(8) :: xi=0d0, log_s=0d0, H_nvt=0d0, temp_var, energy_unit
    real(8), parameter :: T_desired=1d0

    call cpu_time(start_time)
    call RANDOM_SEED(put=seed)
    
    !N_molecules = INT(997.07d3 / 18d0 * N_av * (cell_length * R_0) * (cell_length * R_0) * (cell_length * R_0))
    N_atoms = N_molecules*3
    N_steps = 1d5
    ! N_steps = INT(4d-10/(t_0*dt))

    atom_type(1) = atom(D_null, 1d0, 1d0, -0.82d0)
    atom_type(2) = atom(0.01d0*4.184d3/N_av/E_0, 0.9d-10/R_0, 1d0/16d0, 0.41d0)
    energy_unit = E_0 * E_unit / real(N_atoms, kind=8)

    allocate(type_matrix(N_atoms))
    type_matrix = 2
    do i = 1, N_molecules
        type_matrix(3*(i-1)+1) = 1
    end do
    !print "(3(I0,1X))", type_matrix
    !print *
    allocate(mass(N_atoms), mass_matrix(3, N_atoms), charge(N_atoms))
    mass = atom_type(type_matrix)%mass
    charge = atom_type(type_matrix)%charge
    do i = 1, 3
        mass_matrix(i, :) = mass
    end do
    mass_total = SUM(mass)
    !print "(3(G0,1X))", mass
    !print "(3(G0,1X))", mass_matrix
    !print *, charge

    ! Generating pairs of interaction
    call pairing(N_molecules, pair_bonding, pair_angle, pair_remaining, molecules_comp)
    number_bonds = SIZE(pair_bonding, 2)
    number_angle = SIZE(pair_angle, 2)
    number_interaction = N_atoms * (N_atoms-1) / 2
    number_other = number_interaction - number_bonds - number_angle
    !print "(3(G0,1X))", molecules_comp
    volume = cell_length * cell_length * cell_length
    alpha = (3.6d0 * pi * pi * pi * N_atoms / (volume * volume)) ** (1d0/3d0)       ! For Ewald summation
    rho_null(1) = real(N_molecules, kind=8) / (volume)
    rho_null(2) = real(N_molecules*2, kind=8) / (volume)
    r_cut_third = 1d0 / (r_cut*r_cut*r_cut)
    U_correction = 0d0
    pressure_correction = 0d0
    do i = 1, 2
        do j = 1, 2
            if (i == j) then
                D0_temp = atom_type(i)%D0
                R0_temp = atom_type(i)%R0
            else
                D0_temp = SQRT(atom_type(i)%D0 * atom_type(j)%D0)
                R0_temp = (atom_type(i)%R0 + atom_type(j)%R0) / 2d0
            end if
            U_correction(1) = U_correction(1) + 4d0/3d0 * pi * D0_temp * (R0_temp*R0_temp*R0_temp) * (rho_null(i)*volume) * rho_null(j) * (1d0/6d0 * (r_cut_third*r_cut_third*r_cut_third) - r_cut_third)
            pressure_correction = pressure_correction + 8d0/3d0 * pi * D0_temp * (R0_temp*R0_temp*R0_temp) * rho_null(i) * rho_null(j) * (1d0/3d0 * (r_cut_third*r_cut_third*r_cut_third) - r_cut_third)
        end do
    end do

    allocate(position(3, N_atoms))
    position=0d0
    allocate(velocity(4, N_atoms))
    velocity=0d0
    U_potential=0d0
    U_kinetic=0d0
    U_total=0d0
    pressure = 0d0
    pressure_k = 0d0
    allocate(force(3, N_atoms))
    force=0d0
    xyz=cell_length*R_0*1d9

    k_radius = 2d0
    call find_points_3d(k_radius, k_matrix, k_counts)
    k_matrix = k_matrix * 2d0 * pi / cell_length
    !print "(3(G0,1X))", k_matrix

    U_correction(2) = SQRT(0.5d0 * r_cut / cell_length) / (alpha * r_cut * r_cut) * EXP(- alpha * r_cut * r_cut)
    U_correction(2) = U_correction(2) + SQRT(k_radius * alpha) / (pi * pi * k_radius * k_radius) * EXP(- (pi * pi * k_radius * k_radius) / (alpha * cell_length * cell_length))
    U_correction(2) = U_correction(2) * electric_energy_null * DOT_PRODUCT(charge, charge)

    !call init_position_square(cell_length, N_atoms, position)
    call init_position_square_H2O(cell_length, N_atoms, position, bond_length_eq_null)
    call init_velocity_H2O(N_atoms, mass, velocity)
    !call calc_potential_period(D_null, N_atoms, cell_length, r_cut, position, U_potential, force, pressure)
    call calc_potential_period_LJ(N_atoms, cell_length, r_cut, type_matrix, atom_type, position, U_potential, force, pressure, pair_remaining, number_other)
    call calc_bond(N_atoms, cell_length, position, U_potential, force, pressure, pair_bonding, number_bonds)
    call calc_angle(N_atoms, cell_length, position, U_potential, force, pressure, pair_angle, number_angle)
    call calc_Coulomb(N_atoms, cell_length, r_cut, charge, position, U_potential, force, pressure, pair_remaining, number_other, alpha, k_matrix, k_counts, pair_bonding, number_bonds, pair_angle, number_angle)
    U_kinetic = SUM(velocity(4,:))
    U_total = U_kinetic + U_potential
    T_mean = U_kinetic / real(N_atoms, kind=8)
    pressure = (pressure + U_kinetic*2d0) / (3d0*volume)
    center_velocity = SUM(mass_matrix * velocity(1:3, :), 2) / mass_total

    open(10, file="traj_298K_500.gro", status="Replace")
    open(20, FILE="recording_298K_5000.txt", status="Replace")
    write(10, "('MD of ',I0,' H2O, ','step=', I0, ', time=', G0, ' (ps)')") N_molecules, 0, 0d0
    write(10, "(I0)") N_atoms
    do counter=1, N_molecules
        do i=1, 3
            write(10, '(I5, A5, A5, I5, 3(F8.3), 3(F8.4))') counter, residue_name, atom_name(i), molecules_comp(i, counter), position(1:3, molecules_comp(i, counter))*R_0*1d9, velocity(1:3, molecules_comp(i, counter))*R_0/t_0*1d-3
        end do
    end do
    write(10, '(3F10.5)') xyz

    !write(*,*) "Time scale", t_0
    !write(*,*) "Initial condition"
    !write(*,"(3(F8.3,2X))") position
    !write (*,"(4(F8.4,2X))") velocity
    !write(*,*) "Kinetic energy ", U_kinetic, "Mean: ",U_kinetic/real(N_atoms, kind=8)
    !write(*,*) "Potential energy ", U_potential
    !write(*,*) "Center velocity ", center_velocity
    !print *, "Electric_energy_null = ", electric_energy_null
    !print *, Q_nh_null
    !print *, alpha
    !print *, "sigma = ", SQRT(0.5d0/alpha)

    ! Propagation of position and velocity
    do counter = 1, N_steps
        position = position + velocity*dt + (force / mass_matrix - xi * velocity(1:3,:)) / 2d0 * dt * dt
        position = position - FLOOR(position/cell_length, 8) * cell_length
        velocity(1:3,:) = velocity(1:3,:) + (force / mass_matrix - xi * velocity(1:3,:)) / 2d0 * dt
        temp_var = real(N_atoms, kind=8) / Q_nh_null * (T_mean-T_desired)
        log_s = log_s + (xi + temp_var*dt) * dt
        xi = xi + 2d0 * temp_var * dt
        call calc_potential_period_LJ(N_atoms, cell_length, r_cut, type_matrix, atom_type, position, U_potential, force, pressure, pair_remaining, number_other)
        call calc_bond(N_atoms, cell_length, position, U_potential, force, pressure, pair_bonding, number_bonds)
        call calc_angle(N_atoms, cell_length, position, U_potential, force, pressure, pair_angle, number_angle)
        call calc_Coulomb(N_atoms, cell_length, r_cut, charge, position, U_potential, force, pressure, pair_remaining, number_other, alpha, k_matrix, k_counts, pair_bonding, number_bonds, pair_angle, number_angle)
    
        velocity(1:3,:) = velocity(1:3,:) + (force / mass_matrix - xi * velocity(1:3,:)) / 2d0 * dt
        velocity(4,:) = 0.5d0 * mass * SUM(velocity(1:3,:) * velocity(1:3,:),1)
        U_kinetic = SUM(velocity(4,:))
        U_total = U_kinetic + U_potential
        T_mean = U_kinetic / real(N_atoms, kind=8)
        center_velocity = SUM(mass_matrix * velocity(1:3, :), 2) / mass_total
        pressure_k = (U_kinetic * 2d0) / (3d0 * volume)
        pressure = pressure / (3d0 * volume) + pressure_k     ! preesure_correction needs to be added.
        H_nvt = U_total + SUM(U_correction) + 0.5d0 * xi * xi * Q_nh_null + 2d0 * real(N_atoms, kind=8) * T_desired * log_s
    
        write(*,*)
        write (*,"(A, I0, A, 6(G0, 1X))") 'Step ', counter, ': ', T_mean, U_kinetic, U_potential, pressure_k, pressure, H_nvt
        write(*,"(3(1X,G0))") center_velocity
    
        ! Output to .gro files for VMD visualization
        if (mod(counter, 20)==0) then
            write(20, "(G0,1X,8(G0,1X))") real(counter, kind=8)*dt*t_0, T_mean*T, (pressure)*p_0, (pressure+pressure_correction)*p_0, &
                U_kinetic*energy_unit, U_potential*energy_unit, (U_total)*energy_unit, (U_total+SUM(U_correction))*energy_unit, &
                H_nvt * energy_unit
            write(10, "('MD of ', I0,' H2O, ', 'step=', I0, ', time=', G0, ' (ps)')") N_molecules, counter, counter*t_0*dt*1d12
            write(10, "(I0)") N_atoms
            do i = 1, N_molecules
                do j = 1, 3
                    write(10, '(I5, A5, A5, I5, 3(F8.3), 3(F8.4))') counter, residue_name, atom_name(j), molecules_comp(j, i), position(1:3, molecules_comp(j, i))*R_0*1d9, velocity(1:3, molecules_comp(j, i))*R_0/t_0*1d-3
                end do
            end do
            write(10, '(3F10.5)') xyz
        end if
    
    end do

    close(10)
    close(20)
    deallocate(position, velocity, force)
    deallocate(mass, mass_matrix, k_matrix, charge)
    deallocate(pair_bonding, pair_angle, pair_remaining, molecules_comp, type_matrix)
    call cpu_time(end_time)
    elapsed_time = end_time - start_time
    print*, 'Computation time: ', elapsed_time, ' seconds'


    end program MD_practice_02

