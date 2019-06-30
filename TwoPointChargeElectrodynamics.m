% This script intends to evaluate the effect of the radiated field caused
% by the acceleration on the final energy of a system.

% Many constants have been suppressed, permittivity, and permeability as
% well as the speed of light are equal to 1. 

% 4*pi factor also out, chargues and masses are not subject to any
% particles charge or mass and hence they may take entire values

clear


%initial_pos_1 = [-500, -0.5, 0];
%initial_pos_2 = [500, 0.5, 0];


%initial_vel_1 = [0.1, 0, 0];
%initial_vel_2 = [-0.1, 0, 0];

% Test 1: "1" is an electron and "2" is a massive and stoped proton
%mass1 = 1;
%mass2 = 1000000;

%charge1 = 0.04;
%charge2 = -0.04;

% This test is mainly for frontal approximation, only a small value can be
% added to Y (not add any to Z).
%initial_pos_1 = [-1000, 1, 0];
%initial_pos_2 = [0,0,0];

%initial_vel_1 = [0.2, 0, 0];
%initial_vel_2 = [0, 0, 0];


%Test2: two electrons moving to a frontal colision (that will never happen)

% Particle masses
mass1 = 1;
mass2 = 1;

% Particle charge, same 
charge1 = 1;
charge2 = 1;    

% Large distance, the bigger the distance is, the bigger are the intervals
% evaluated, distance intervals are recalculated based on the distance
initial_pos_1 = [-500000000, 0, 0];
initial_pos_2 = [500000000, 0, 0];

% Particle velicity, it should be on X axis, the speed can be different 
% for each particle
initial_vel_1 = [0.1, 0, 0];
initial_vel_2 = [-0.1, 0, 0];


% prebuilt variables to reduce the amount of calculus.
charge_product = charge1 * charge2;
charge1_square = charge1^2;
charge2_square = charge2^2;


% It approximates how much the particles may get closer or farther
difference_velocity = norm(initial_vel_1 - initial_vel_2 );

initial_distance = norm(initial_pos_1 - initial_pos_2);

% Number of iterations before the initial gets multiplied or divided by 2
num_it_segment = 2000;                % Max used: 24000 

dr_iteration = initial_distance / num_it_segment / 2;

iteration_dif_time = dr_iteration / difference_velocity;

% Forecast the vector size 
estimated_segment_number = round(log2(norm(initial_pos_2 - initial_pos_1))) + 10;      


%vector_size = 2^(estimated_segment_number) * num_it_segment  * 2;
vector_size = estimated_segment_number * num_it_segment * 2;


position_1 = zeros(vector_size, 3);
position_2 = zeros(vector_size, 3);
velocity_1 = zeros(vector_size, 3);
velocity_2 = zeros(vector_size, 3);
% Calculate also the momentum
mommentum_1 = zeros(vector_size, 3);
mommentum_2 = zeros(vector_size, 3);
iteration_time = zeros(vector_size, 3);


acceleration_1 = zeros(vector_size, 3);
acceleration_2 = zeros(vector_size, 3);
larmor_energy_1 = 0;
larmor_energy_2 = 0;


position_1(1,:) = initial_pos_1;
position_2(1,:) = initial_pos_2;
velocity_1(1,:) = initial_vel_1;
velocity_2(1,:) = initial_vel_2;
mommentum_1(1,:) = mass1 / sqrt( 1 - norm(velocity_1(1,:))^2) * velocity_1(1,:);
mommentum_2(1,:) = mass2 / sqrt( 1 - norm(velocity_2(1,:))^2) * velocity_2(1,:);


iter = 2;       
iter_old = 1;   
time = 0;
% This is done to force the retarded position to be calculated 
% assuming constant velocity
retarded_time_1 = -1;
retarded_time_2 = -1;

distance = initial_distance;
original_initial_distance = initial_distance;      % initial_distance is updated in each iteration!

maximmum_velocity_1 = 0;
maximmum_velocity_2 = 0;
minimmum_velocity_1 = 1;
minimmum_velocity_2 = 1;
speed_1 = norm(velocity_1(1,:));
speed_2 = norm(velocity_2(1,:));

initial_potential_1 = 0;

while iter < vector_size && distance < (original_initial_distance + 1)
    time = time + iteration_dif_time;
    
    position_1(iter,:) = position_1(iter_old,:) + velocity_1(iter_old,:) * iteration_dif_time;
    position_2(iter,:) = position_2(iter_old,:) + velocity_2(iter_old,:) * iteration_dif_time;
  
    
% Obtain the retarded distance ignoring Y abd Z because the effect will be
% minimal.
    if retarded_time_1 < 0
        retarded_time_1 = time - distance / (1- norm(velocity_2(1,:)));
        retarded_position_1 = position_2(iter,:) - velocity_2(iter_old,:) * (time - retarded_time_1);
        
        if retarded_time_1 >= 0
            retarded_position_index_1 = 1;
        end
        retarded_velocity_1 = velocity_2(1,:);
        
        retarded_acceleration_1 = [0,0,0];  % Ignore the retarded acceleration when far away
     
% This potential is Potencial will only be calculated on "1", for
% considereing this potencial, use same speed on both particles
        if initial_potential_1 == 0
            [initial_potential_1, initial_vector_potential_1] = Lienard_Wiechert_Potential( (position_1(iter,:)-retarded_position_1), initial_vel_2);
            initial_potential_1 = initial_potential_1 * charge_product;
            initial_vector_potential_1 = initial_vector_potential_1 * charge_product;
            initial_lagrangian_1 = initial_potential_1 - dot(initial_vector_potential_1, initial_vel_1);
        end
        
% Algorithm to get the retarded position 
    else
       do_continue = 1;
       while do_continue == 1
          retarded_position_index_1_next = retarded_position_index_1 + 1;
          auxiliary_distance_1 = norm( position_1(iter,:) - position_2(retarded_position_index_1_next,:)) ;
          
% Check if the ligth from the previous position is later than the next
% retarded field position, if this is not the case we should move on.
          if ( time - iteration_time(retarded_position_index_1_next) > auxiliary_distance_1 )
            retarded_position_index_1 = retarded_position_index_1 + 1;  
          else
            do_continue = 0;  
          end        
       end
       interval = iteration_time(retarded_position_index_1_next) - iteration_time(retarded_position_index_1);
       
       cof = abs(time - iteration_time(retarded_position_index_1) - auxiliary_distance_1) / interval;
       
% Do a lineal interpolation on the initial position and velocity 
       retarded_position_1 = (1-cof) * position_2(retarded_position_index_1,:) + cof * position_2(retarded_position_index_1_next,:); 
       retarded_velocity_1 = (1-cof) * velocity_2(retarded_position_index_1,:) + cof * velocity_2(retarded_position_index_1_next,:);
       retarded_acceleration_1 = (1-cof) * acceleration_2(retarded_position_index_1,:) + cof * acceleration_2(retarded_position_index_1_next,:);
       
    end
    
    
    if retarded_time_2 < 0
        retarded_time_2 = time - distance / (1- norm(velocity_1(1,:)));
        retarded_position_2 = position_1(iter,:) - velocity_1(iter_old,:) * (time - retarded_time_2);
        
        if retarded_time_2 >= 0
            retarded_position_index_2 = 1;
        end
        retarded_velocity_2 = velocity_1(1,:);
        
        retarded_acceleration_2 = [0,0,0];  % Ignore the retarded acceleration when far away
        
% Algorithm to get the retarded position    
    else
       do_continue = 1;
       while do_continue == 1
          retarded_position_index_2_next = retarded_position_index_2 + 1;
          auxiliary_distance_2 = norm( position_2(iter,:) - position_1(retarded_position_index_2_next,:)) ;
          
% Check if the ligth from the previous position is later than the next
% retarded field position, if this is not the case we should move on.
          if ( time - iteration_time(retarded_position_index_2_next) > auxiliary_distance_2 )
            retarded_position_index_2 = retarded_position_index_2 + 1;  
          else
            do_continue = 0;  
          end        
       end
       interval = iteration_time(retarded_position_index_2_next) - iteration_time(retarded_position_index_2);
       cof = abs(time - iteration_time(retarded_position_index_2) - auxiliary_distance_2) / interval;
       
% Do a lineal interpolation on the initial position and velocity 
       retarded_position_2 = (1-cof) * position_1(retarded_position_index_2,:) + cof * position_1(retarded_position_index_2_next,:);
       retarded_velocity_2 = (1-cof) * velocity_1(retarded_position_index_2,:) + cof * velocity_1(retarded_position_index_2_next,:);
       retarded_acceleration_2 = (1-cof) * acceleration_1(retarded_position_index_2,:) + cof * acceleration_1(retarded_position_index_2_next,:);
       
    end    
    
% Get EM field and Lorentz force 
    retarded_vector_1 = position_1(iter,:) - retarded_position_1;
    retarded_vector_2 = position_2(iter,:) - retarded_position_2;
    
    [E1, B1] = Induced_Field_Without_Units(retarded_vector_1, retarded_velocity_1);
    [E2, B2] = Induced_Field_Without_Units(retarded_vector_2, retarded_velocity_2);    
    
    [E_rad_1, B_rad_1] = RadiatedFieldAcceleratedCharge( retarded_vector_1, retarded_velocity_1, retarded_acceleration_1);
    [E_rad_2, B_rad_2] = RadiatedFieldAcceleratedCharge( retarded_vector_2, retarded_velocity_2, retarded_acceleration_2);
     
    E1 = E1 + E_rad_1;
    B1 = B1 + B_rad_1;
    E2 = E2 + E_rad_2;
    B2 = B2 + B_rad_2;  

    
    F1 = Lorentz_Force( velocity_1(iter,:), E1, B1);
    F2 = Lorentz_Force( velocity_2(iter,:), E2, B2);   
    
   
    dif_momentum_1 = F1 * iteration_dif_time * charge_product;
    dif_momentum_2 = F2 * iteration_dif_time * charge_product;
   
    sigma_1 = 1 / sqrt( 1 - speed_1^2);
    sigma_2 = 1 / sqrt( 1 - speed_2^2);
        
% Acceleration measured in the lab, it is a rough approximation but
% will be good because velocity is post on  X direction.
    acceleration_1(iter,:) = [dif_momentum_1(1) / sigma_1, dif_momentum_1(2), dif_momentum_1(3)] / (mass1 * sigma_1 * iteration_dif_time);
    acceleration_2(iter,:) = [dif_momentum_2(1) / sigma_2, dif_momentum_2(2), dif_momentum_2(3)] / (mass2 * sigma_2 * iteration_dif_time);
      
    larmor_energy_1 = larmor_energy_1 + norm(acceleration_1(iter,:))^2 * (2/3) * charge1_square * iteration_dif_time;
    larmor_energy_2 = larmor_energy_2 + norm(acceleration_1(iter,:))^2 * (2/3) * charge2_square * iteration_dif_time;
    
    mommentum_1(iter,:) = mommentum_1(iter_old,:) + dif_momentum_1;
    mommentum_2(iter,:) = mommentum_2(iter_old,:) + dif_momentum_2;
        
% v^2(1+p^2) = p^2 / m^2   
    mom_1_cuadr = norm(mommentum_1(iter,:))^2;
    mom_2_cuadr = norm(mommentum_2(iter,:))^2;
    v_1_abs = sqrt( mom_1_cuadr / (1 + mom_1_cuadr) ) / mass1;
    v_2_abs = sqrt( mom_2_cuadr / (1 + mom_2_cuadr) ) / mass2;
    
    velocity_1(iter,:) = v_1_abs * mommentum_1(iter,:) / norm(mommentum_1(iter,:));
    velocity_2(iter,:) = v_2_abs * mommentum_2(iter,:) / norm(mommentum_2(iter,:)); 
    
    speed_1 = norm(velocity_1(iter,:));
    speed_2 = norm(velocity_2(iter,:));
    
    if speed_1 > maximmum_velocity_1 
        maximmum_velocity_1 = speed_1;
    elseif speed_1 < minimmum_velocity_1 
        minimmum_velocity_1 = speed_1;
    end
    
    if speed_2 > maximmum_velocity_2 
        maximmum_velocity_2 = speed_2;
    elseif speed_2 < minimmum_velocity_2 
        minimmum_velocity_2 = speed_2;
    end    
    
% Last, check if a change in resolution is needed.
    current_distance = norm(position_1(iter,:) - position_2(iter,:));
    
    if current_distance <= (initial_distance / 2) || current_distance >= (initial_distance * 2)
        initial_distance = current_distance;
        dr_iteration = initial_distance / num_it_segment / 2;
        difference_velocity = norm(velocity_1(iter,:) - velocity_2(iter,:));
        iteration_dif_time = dr_iteration / difference_velocity;  
    end
    
    
    iter_old = iter;
    iteration_time(iter) = time; 
    distance = norm( position_1(iter,:) - position_2(iter,:));
    iter = iter + 1;
    
end


% Show results
numer_of_iterations = iter_old

final_position_1= position_1(iter_old,:)
final_position_2= position_2(iter_old,:)
final_velocity_1 = velocity_1(iter_old,:)
final_velocity_2 = velocity_2(iter_old,:)

%speed_difference_1 = norm(final_velocity_1) - norm(initial_vel_1)
%speed_difference_2  = norm(final_velocity_2) - norm(initial_vel_2)

final_sigma_1 = 1 / sqrt(1 - norm(final_velocity_1)^2);
final_sigma_2 = 1 / sqrt(1 - norm(final_velocity_2)^2);

initial_sigma_1 =  1 /sqrt(1 - norm(initial_vel_1)^2);
initial_sigma_2 =  1 /sqrt(1 - norm(initial_vel_2)^2);


kinetic_energy_difference_1 = mass1 * (final_sigma_1 - initial_sigma_1) 
kinetic_energy_difference_2 = mass2 * (final_sigma_2  - initial_sigma_2)


maximmum_velocity_1 = maximmum_velocity_1
maximmum_velocity_2 = maximmum_velocity_2
minimmum_velocity_1 = minimmum_velocity_1
minimmum_velocity_2 = minimmum_velocity_2


% This is only valid when both particles have same mass and charge and
% opossite initial velocity
minimmum_sigma_1 = 1 / sqrt(1 - minimmum_velocity_1^2);
minimmum_energy_1 = (minimmum_sigma_1 - 1) * mass1;

initial_energy_1 = (initial_sigma_1 - 1) * mass1;

[V_final_1, A_final_1] = Lienard_Wiechert_Potential( (position_1(iter_old,:)-retarded_position_1), retarded_velocity_1);
V_final_1 = V_final_1 * charge_product
A_final_1 = A_final_1 * charge_product;

difference_of_potential = V_final_1 - initial_potential_1

%Lagr_final_1 = V_final_1 - dot(A_final_1, final_velocity_1)
%Dif_lagr = Lagr_fin_1 - initial_lagrangian_1

% Show the accumulated larmor energy.
larmor_energy_1 = larmor_energy_1 
larmor_energy_2 = larmor_energy_2

