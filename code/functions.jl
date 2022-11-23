""" 
Create an hexagonal lattice by placing a node at each vertex.
    nr = number of particles in row
    nt = number fo particles in width
    dh = horizontal space between particles
    dv = verticle space between particles
    reg = if true the hexagons will be regular, if false we introduce a perturbation to the original regular structure
    pert = amount of perturbation to add
"""
function PositionParticles(nr::Integer, nt::Integer, dh::Real, dv::Real; reg::Bool=true, pert = 0.5 )

    pos = Array{Float64}(undef,0,2)
    x =  collect(0.0:dh:(dh*(nr-1)))
    h = 0;
    id = 1
    for i=1:1:nt
      if isodd(i)
        for j = 1:1:nr
          new = [x[j] h]
          pos = vcat( pos, new)
          id = id+1
        end
      else
        for j = 1:1:nr-1
          new = [x[j].+dh/2 h]
          pos = vcat( pos, new)
          id = id + 1
        end
      end
      h = h+dv;
    end
  
    if reg ==false
    # If we want a random geometry: we introduce a perturbation to the original regular structure.
      for i =1:1:size(pos,1)
        pos[i,:] = pos[i,:] .+ rand(2)*dv*pert
      end
    end

    # Vector containing the qs relative to each node
    nodes_q = zeros(size(pos,1),2)
    count = 1;
    for i = 1:1:size(pos,1)
        nodes_q[i,1] = count;
        nodes_q[i,2] = count+1;
        count = count + 2;
    end
    nodes_q = convert(Array{Int64},nodes_q)

    connectivity = LinkParticles(nr,nt)


    return pos, nodes_q, connectivity
end



""" 
Create the connections between particles.    
Function called by PositionParticles()
"""
function LinkParticles(nr, nt)
    if iseven(nt)
      n_link = (nt/2) * (2*nr-3) + (nr-1)*2 *(nt-1)
    else
      n_link = ceil(nt/2) * (nr-1) + (ceil(nt/2)-1) * (nr-2) + (nr-1)*2*(nt-1)
    end
    connectivity = Array{Float64}(undef, 0, 2)
    # ---------------------Create horizontal links------------------------------
    index = 1;
    for i = 1:1:nt
      isodd(i) ? max = nr-1 : max = nr-2
      for j = 1:1:max
          new = [index index+1]
          connectivity = vcat(connectivity, new)
          index = index + 1;
      end
      index = index +1;
     end
    # ---------------------Create tilted links----------------------------------
    index = 1;
    for i = 1:1:(nt-1)
      if isodd(i)
        index = index + nr;
        for j = 1:1:(nr-1)
            new = [index-nr index]
            connectivity = vcat(connectivity, new)
            new = [index-(nr-1) index]
            connectivity = vcat(connectivity, new)
          index = index+1;
        end
      elseif iseven(i)
        index = index- (nr-1);
        for j = 1:1:(nr-1)
            new = [index index+(nr-1)]
            connectivity = vcat(connectivity, new)
            new = [index index+nr]
            connectivity = vcat(connectivity, new)
          index = index+1;
        end
      end
    end

    connectivity = convert(Array{Int64},connectivity)

    return connectivity
end


""" 
Select the cells within an area identified by x1, x2, y1, y2.    
"""
function SelectCell(position, connectivity, x1, x2, y1, y2)

    n_nodes = Vector{Integer}(undef,0)
    e = 1e-4
        fx_index = findall(x->(x2+e)>=x && (x1-e)<=x ,position[:,1])
        fy_index = findall(y->(y1-e)<=y && (y2+e)>=y,position[:,2])
        final_index = findall(in(fy_index), fx_index)
        indexes = fx_index[final_index]
  
        e_cell = Vector{Integer}(undef,0)
        for i = 1:1:size(indexes,1)
  
          all_e1 = findall(x->x==indexes[i] , connectivity[:,:])
          all_e1 = map(i->i[1], all_e1)
  
          for j = 1:1:size(indexes,1)
  
            if j != i
               all_e2 = findall(x->x==indexes[j] , connectivity[:,:])
               all_e2 = map(i->i[1], all_e2)
               final_e = findall(in(all_e1), all_e2)
               e_cell = vcat(e_cell,all_e2[final_e])
            end
  
          end
  
        end
  
    return indexes, e_cell
  
end

"""
 Function to generate the global sfiffness matrix of the system
 it requires the positions (x,y) of the nodes; the relative qs for each node (e.g. q1,q2),
 the connectivity (1,2, where 1 2 are the nodes numbers) and the stiffness of the bar. 
"""
function K_global(nodes_position, nodes_q, connectivity, EA, prestress; selected = [], EAmultiplier = nothing, preMultiplier = nothing)

    # Find angles rotation from local coordinates to global coordinates and lenght
    geo = findAngles(nodes_position, connectivity);

    n = size(nodes_position,1);
    A = zeros(n*2,n*2);
    F = zeros(n*2);

    for i = 1:1:size(connectivity,1)

        found = findall(x->x==i, selected)
        k = EA
        pre = prestress
        if !isempty(found)
          k = EA *EAmultiplier;
          pre = prestress*preMultiplier;
        end

        q1 = nodes_q[connectivity[i,1],1];
        q2 = nodes_q[connectivity[i,1],2];

        q3 = nodes_q[connectivity[i,2],1];
        q4 = nodes_q[connectivity[i,2],2];

        # Use for debugging
                # print("i: ", i, "\n");
                # print("connectivity: ", connectivity[i,:], "\n");

                # print("q1: ", q1, " ");
                # print("q2: ", q2, " ");
                # print("q3: ", q3, " ");
                # print("q4: ", q4, "\n");

        A = assembleMatrix(A, k, geo[i,:], q1, q2, q3, q4);
        F = assembleForce(F, pre, geo[i,:], q1, q2, q3, q4)

    end

    return A,F
end

"""
 Function to calculate the cos, sin and length of each element. needed to apply the rotation to the 
 element stiffness matrix
 Parameters needed: positions of the nodes (x,y) and the connectivity of each element
"""
function findAngles(nodes, connectivity)
    c = zeros(size(connectivity,1))
    s = zeros(size(connectivity,1))
    l = zeros(size(connectivity,1))

    for i = 1:1:size(c,1)

        p1 = connectivity[i,1];
        p2 = connectivity[i,2];
        x1 = nodes[p1,1];
        y1 = nodes[p1,2];
        x2 = nodes[p2,1];
        y2 = nodes[p2,2];
        L = sqrt((x2-x1)^2 + (y2-y1)^2);

        c[i] = (x2-x1)/L;
        s[i] = (y2-y1)/L;
        l[i] = L;

    end

    geo = [c s l]

    return geo

end

"""
 Function to assemble the element stiffness matrix into the global one.
 Parameters needed: the globa stiffness matrix, the stiffness of the bar (EA), the qs that the element connects
"""
function assembleMatrix(A, k, geo, q1, q2, q3, q4)

    # Element stiffness matrix
    c = geo[1];
    s = geo[2];
    l = geo[3];

    a = (k/l)* [c^2   c*s  -c^2 -c*s
                c*s   s^2  -c*s -s^2
                -c^2  -c*s  c^2  c*s
                -c*s  -s^2  c*s  s^2]

    # Pasting the element stiffness matrix into the Global stiffness matrix
    A[q1,q1] = A[q1,q1] + a[1,1];
    A[q1,q2] = A[q1,q2] + a[1,2];
    A[q1,q3] = A[q1,q3] + a[1,3];
    A[q1,q4] = A[q1,q4] + a[1,4];     

    A[q2,q1] = A[q2,q1] + a[2,1];
    A[q2,q2] = A[q2,q2] + a[2,2];
    A[q2,q3] = A[q2,q3] + a[2,3];
    A[q2,q4] = A[q2,q4] + a[2,4]; 

    A[q3,q1] = A[q3,q1] + a[3,1];
    A[q3,q2] = A[q3,q2] + a[3,2];
    A[q3,q3] = A[q3,q3] + a[3,3];
    A[q3,q4] = A[q3,q4] + a[3,4]; 

    A[q4,q1] = A[q4,q1] + a[4,1];
    A[q4,q2] = A[q4,q2] + a[4,2];
    A[q4,q3] = A[q4,q3] + a[4,3];
    A[q4,q4] = A[q4,q4] + a[4,4]; 

    return A

end

"""
Function to assemble the forces at the nodes.
"""
function assembleForce(F, prestress, geo, q1, q2, q3, q4)

    # Element stiffness matrix
    c = geo[1];
    s = geo[2];
    l = geo[3];
  
    f = prestress* [c -s
                     s  c ]
  
    # Pasting the element stiffness matrix into the Global stiffness matrix
    F[q1] = F[q1] + f[1];    
  
    F[q2] = F[q2] + f[2];
  
    F[q3] = F[q3] - f[1];
  
    F[q4] = F[q4] - f[2];
  
    return F
  
end


"""
Select all nodes with a given x location or y location
"""
function SelectNodes(position; x=NaN,y=NaN)
    n_nodes = Vector{Integer}(undef,0)
    e = 1e-3
      if !isnan(x) && isnan(y)
        f_index = findall(xx->(xx+e)>=x && (xx-e)<=x ,position[:,1])
      elseif isnan(x) && !isnan(y)
        f_index = findall(yy->(yy+e)<=y && (yy-e)>=y,position[:,2])
      end
    return f_index
end

"""
 Delete rows and columns of the contrained nodes (node_to_fix) from the global stiffness matrix (K), 
 the Force and the q vector (it will be important to know what qs are the one calculated to update the relative values). 
 Displacement BC are applied by removing from the force vector the [column from K of the q] * disp value. 
"""
function applyBC(K, F_vect, node_right, node_left, nodes_q, disp; middle = nothing, middle_length =nothing)
  
    q_vect_disp_right = Array{Int64}(undef,0);
    for i = 1:1:size(node_right,1)
        push!(q_vect_disp_right,nodes_q[node_right[i],1])
    end

    for i = 1:1:size(q_vect_disp_right,1)

      F_vect = F_vect - disp * K[:,q_vect_disp_right[i]];

    end


    q_vect_disp_left = Array{Int64}(undef,0);
    for i = 1:1:size(node_left,1)
        push!(q_vect_disp_left,nodes_q[node_left[i],1])
    end

    for i = 1:1:size(q_vect_disp_left,1)

      F_vect = F_vect - (-disp) * K[:,q_vect_disp_left[i]];

    end


    q_vect_fix = Array{Int64}(undef,0);
    for i = 1:1:size(node_right,1)
      push!(q_vect_fix,nodes_q[node_right[i],2])
    end

    for i = 1:1:size(node_left,1)
      push!(q_vect_fix,nodes_q[node_left[i],2])
    end

    if !isnothing(middle)

      for i = 1:1:size(node_middle,1)
        push!(q_vect_fix,nodes_q[node_middle[i],1])
      end
    end

    if !isnothing(middle_length)

      for i = 1:1:size(node_middle_length,1)
        push!(q_vect_fix,nodes_q[node_middle_length[i],2])
      end

    end


    q_vect = sort([q_vect_disp_right; q_vect_disp_left; q_vect_fix])
    q_fix = Tuple(Int64(x) for x in q_vect)

    q_to_update = collect(1:size(K,1));

    newK = K[setdiff(1:end, q_fix), setdiff(1:end, q_fix)];
    newF = F_vect[setdiff(1:end, q_fix)];
    newq_to_update = q_to_update[setdiff(1:end, q_fix)];
    return newK, newF, newq_to_update, q_vect_disp_right, q_vect_disp_left

end

"""
 Given the new displacements we can caluclate the new node positions
 Parameters needed: old node positions, the displacement vector (q) and the vectors of the qs that have been calculated.
"""
function update_position(old_position, nodes_q, q, q_update, q_disp, q_disp_left, disp)

    new_position = copy(old_position)

    for i = 1:1:size(q,1)

        node_q1 = findall(x -> x == q_update[i], nodes_q[:,1])
        node_q2 = findall(x -> x == q_update[i], nodes_q[:,2])

        if !isempty(node_q1)
            new_position[node_q1[1],1] = new_position[node_q1[1],1] + q[i];
        elseif !isempty(node_q2)
            new_position[node_q2[1],2] = new_position[node_q2[1],2] + q[i];
        end

    end

    for i = 1:1:size(q_disp,1)

      node_q1 = findall(x -> x == q_disp[i], nodes_q[:,1])
      node_q2 = findall(x -> x == q_disp[i], nodes_q[:,2])

      if !isempty(node_q1)
          new_position[node_q1[1],1] = new_position[node_q1[1],1] + disp;
      # elseif !isempty(node_q2)
      #     new_position[node_q2[1],2] = new_position[node_q2[1],2] + disp;
      end

    end

    for i = 1:1:size(q_disp_left,1)

      node_q1 = findall(x -> x == q_disp_left[i], nodes_q[:,1])
      node_q2 = findall(x -> x == q_disp_left[i], nodes_q[:,2])

      if !isempty(node_q1)
          new_position[node_q1[1],1] = new_position[node_q1[1],1] - disp;
      end

    end


    all_q = Array{Int64}(undef,0);

    for i = 1:1:size(old_position,1)

      all_q = vcat(all_q, new_position[i,1]-old_position[i,1]);
      all_q = vcat(all_q, new_position[i,2]-old_position[i,2]);

    end

    return new_position, all_q

end






function calculate_tensor(selectedcell, nodes_position, connectivity, stress)


  elements01 = findall(x->x==selectedcell[1], connectivity[:,1])
  elements02 = findall(x->x==selectedcell[1], connectivity[:,2])
  elements = sort(vcat(elements01, elements02))

  q_hexagon = Array{Int64}(undef,0,1)
  for i = 1:1:size(elements01,1)
    q_hexagon = vcat(q_hexagon, connectivity[elements01[i],2])

  end
  for i = 1:1:size(elements02,1)
    q_hexagon = vcat(q_hexagon, connectivity[elements02[i],1])

  end

  coordinates_q_hexagon = Array{Float64}(undef,0,2)
  for i = 1:1:size(elements,1)
    coordinates_q_hexagon = vcat(coordinates_q_hexagon,nodes_position[q_hexagon[i],:]')
  end
  

  center_coordinates = nodes_position[selectedcell,:];
  shift_coordinates = copy(coordinates_q_hexagon);
  shift_coordinates[:,1] = shift_coordinates[:,1] .- center_coordinates[1,1];
  shift_coordinates[:,2] = shift_coordinates[:,2] .- center_coordinates[1,2];
  ordered_coordinates = shift_coordinates[1,:]';

  #print("\n shifted coordinates", shift_coordinates, "\n")

  angles = zeros(size(shift_coordinates,1),2);

  for i = 1:1:size(shift_coordinates,1)

      vect_mod = sqrt(shift_coordinates[i,1]^2 + shift_coordinates[i,2]^2);

      cosa = shift_coordinates[i,1]/vect_mod;

      angles[i,1] = atan(shift_coordinates[i,2],shift_coordinates[i,1]);
      angles[i,2] = i;
  end
  angles = sortslices(angles, dims = 1, by = x -> x[1])

  ordered_coordinates = zeros(size(shift_coordinates,1),2)

  for i = 1:1:size(angles,1)

      ordered_coordinates[i,:] = shift_coordinates[convert(Int64,angles[i,2]),1:2];
  end

  ordered_coordinates[:,1] = ordered_coordinates[:,1] .+ center_coordinates[1,1];
  ordered_coordinates[:,2] = ordered_coordinates[:,2] .+ center_coordinates[1,2];

  clockwise_sum = 0

  for i = 1:1:size(ordered_coordinates,1)-1

    clockwise_sum = clockwise_sum + (ordered_coordinates[i,1]-ordered_coordinates[i+1,1]) * (ordered_coordinates[i,2]-ordered_coordinates[i+1,2]);

  end

  if clockwise_sum < 0

    ordered_coordinates = reverse(ordered_coordinates,dims = 1)
  end



  n = size(shift_coordinates,1)
  area = 0
  for i = 1:1:n-1

    area = area + ordered_coordinates[i,1] * ordered_coordinates[i+1,2] - ordered_coordinates[i,2] * ordered_coordinates[i+1,1]
    
  end

  area = area + ordered_coordinates[end,1] * ordered_coordinates[1,2] - ordered_coordinates[end,2] * ordered_coordinates[1,1]

  area = abs(area)

  #plot(ordered_coordinates[:,1], ordered_coordinates[:,2])

  geo = findAngles(nodes_position, connectivity);

  tensor = zeros(2,2)
  l_elements = [];

  for i = 1:1:size(elements,1)
    n = elements[i];

    c = geo[n,1];
    s = geo[n,2];
    l = geo[n,3]./2;
    lx = l*c;
    ly = l*s;
    l_elements = vcat(l_elements,l)

    tensor = tensor .+ stress[n]./l .* [lx; ly] * [lx ly];

  end

  min_x = minimum(l_elements);
  max_x = maximum(l_elements);

  A = pi*min_x *max_x
  print("Area = ", area/4, "\n");

  tensor = tensor ./ (area/4);

  display("text/plain", tensor)
  print("\n")

  return(tensor)

end






function calculate_strain(selectedcell, nodes_position, connectivity, strain)


  elements01 = findall(x->x==selectedcell[1], connectivity[:,1])
  elements02 = findall(x->x==selectedcell[1], connectivity[:,2])
  elements = sort(vcat(elements01, elements02))

  for i = 1:1:size(elements,1)

    n = elements[i];
    print("Strain = ", strain[n], "\n")


  end


end






# -------------------------------------------------------------
#                  Plotting functions
#--------------------------------------------------------------

""" 
Plot lattice geometry    
"""
function plot_lattice(node_positions, connectivity, color, linestyle; selected = nothing)

    (fig2, ax2) = subplots(1, 1, figsize=(7,5))
     ax2.set_facecolor((0, 0, 0))
    # ax2.plot(node_positions[:,1], node_positions[:,2], "o", color = color);
 
     for i = 1:1:size(connectivity,1)
 
         x1 = node_positions[connectivity[i,1],1]
         y1 = node_positions[connectivity[i,1],2]
 
         x2 = node_positions[connectivity[i,2],1]
         y2 = node_positions[connectivity[i,2],2]
 
 
         ax2.plot([x1,x2], [y1,y2], color = color, linestyle = linestyle);
 
     end
 
 
     if !isnothing(selected)
 
         for i = 1:1:size(selected,1)
 
             x1 = node_positions[connectivity[selected[i],1],1]
             y1 = node_positions[connectivity[selected[i],1],2]
 
             x2 = node_positions[connectivity[selected[i],2],1]
             y2 = node_positions[connectivity[selected[i],2],2]
 
 
             ax2.plot([x1,x2], [y1,y2], color = "red", linestyle = linestyle);
 
         end
     end
 
 
end




function plot_stress(node_positions, new_positions, EA, connectivity, prestress; selected = nothing, EAmultiplier = nothing, preMultiplier = nothing, stress_noinclusion = nothing)

  σ = zeros(size(connectivity,1))

   for i = 1:1:size(connectivity,1)

       x1 = node_positions[connectivity[i,1],1]
       y1 = node_positions[connectivity[i,1],2]

       x2 = node_positions[connectivity[i,2],1]
       y2 = node_positions[connectivity[i,2],2]

       L0 = sqrt((x2-x1)^2+(y2-y1)^2);

       x1 = new_positions[connectivity[i,1],1]
       y1 = new_positions[connectivity[i,1],2]

       x2 = new_positions[connectivity[i,2],1]
       y2 = new_positions[connectivity[i,2],2]

       L1 = sqrt((x2-x1)^2+(y2-y1)^2);

       σ[i] = EA/L0 *(L1-L0) + prestress;

   end


   if !isnothing(selected)

       for i = 1:1:size(selected,1)

           x1 = node_positions[connectivity[selected[i],1],1]
           y1 = node_positions[connectivity[selected[i],1],2]

           x2 = node_positions[connectivity[selected[i],2],1]
           y2 = node_positions[connectivity[selected[i],2],2]

           L0 = sqrt((x2-x1)^2+(y2-y1)^2);

           x1 = new_positions[connectivity[selected[i],1],1]
           y1 = new_positions[connectivity[selected[i],1],2]
    
           x2 = new_positions[connectivity[selected[i],2],1]
           y2 = new_positions[connectivity[selected[i],2],2]

           L1 = sqrt((x2-x1)^2+(y2-y1)^2);

           σ[selected[i]] = (EA *EAmultiplier)/L0 *(L1-L0) + prestress*preMultiplier;

       end
   end

  
  (fig1, ax1) = subplots(1, 1, figsize=(12,12))
  #(fig1, ax1) = subplots(1, 1, figsize=(7,5)) # Figures papers
  ax1.set_facecolor((0, 0, 0))

  max_σ = maximum(σ);
  min_σ = minimum(σ);




  for i = 1:1:size(connectivity,1)

    x1 = new_positions[connectivity[i,1],1]
    y1 = new_positions[connectivity[i,1],2]

    x2 = new_positions[connectivity[i,2],1]
    y2 = new_positions[connectivity[i,2],2]

    if σ[i] >= max_σ
      σ[i] = max_σ
    elseif σ[i] <= min_σ
      σ[i] = min_σ
    end

    if σ[i] <= 0
      b = 1
      r = -σ[i]/ min_σ + 1
      g = -σ[i]/ min_σ + 1
    elseif σ[i] > 0
      r = 1
      g = -σ[i]/ max_σ + 1
      b = -σ[i]/ max_σ + 1
    end


    ax1.plot([x1,x2], [y1,y2], linestyle = "-", color = (r,g,b), linewidth = 3);
  end
  

  xmax = maximum(new_positions[:,1]) + 0.2
  ymax = maximum(new_positions[:,2]) + 0.2
  xmin = minimum(new_positions[:,1]) - 0.2
  ymin = minimum(new_positions[:,2]) - 0.2

  print("max x = ", xmax, "\n")
  print("max y = ", ymax, "\n")
  xmax = 5.4;
  ymax = 4.4;

  xmin = -1.2;
  ymin = -0.4;

  ax1.set_xlim([xmin, xmax])
  ax1.set_ylim([ymin, ymax])

   ax1.set_xticklabels([])
   ax1.set_xticks([])
   ax1.set_yticklabels([])
   ax1.set_yticks([])

  line_z = collect(min_σ:0.001:max_σ) 
  stepy = ymax/(length(line_z))
  line_y = collect(0:stepy:ymax)

  for i = 1:1:size(line_z,1)-1

    if line_z[i] <= 0
      b = 1
      r = -line_z[i]/ min_σ + 1
      g = -line_z[i]/ min_σ + 1
    elseif line_z[i] > 0
      r = 1
      g = -line_z[i]/ max_σ + 1
      b = -line_z[i]/ max_σ + 1
    end

  
    ax1.plot([xmax, xmax], [line_y[i], line_y[i+1]], linestyle = "-", color = (r,g,b), linewidth = 16);
  end

  annotate(round(min_σ, digits = 3),
            xy=[1;0],
            xycoords="axes fraction",
            xytext=[-30,5],
            textcoords="offset points",
            fontsize=10.0,
            ha="right",
            va="bottom",
            color = "white") 

  annotate(round(max_σ, digits = 3),
          xy=[1;1],
          xycoords="axes fraction",
          xytext=[-30,-5],
          textcoords="offset points",
          fontsize=10.0,
          ha="right",
          va="top",
          color = "white") 
  
  ax1.set_aspect("equal")
  
  
  
  if !isnothing(stress_noinclusion)

    selectedcell, e_selected = SelectCell(node_positions, connectivity, 1.5, 2.5, 1.5, 2.5)


    (fig2, ax2) = subplots(1, 1, figsize=(7,5))
    ax2.set_facecolor((0, 0, 0))
    ratio_σ = σ - stress_noinclusion;
    ratio_plot = []
    for i = 1:1:size(e_selected,1)
      ratio_plot = vcat(ratio_plot, ratio_σ[e_selected[i]])
    end
  
    max_ratio_σ = maximum(ratio_plot);
    min_ratio_σ = minimum(ratio_plot);
  
    for i = 1:1:size(e_selected,1)
  
      x1 = node_positions[connectivity[e_selected[i],1],1]
      y1 = node_positions[connectivity[e_selected[i],1],2]
  
      x2 = node_positions[connectivity[e_selected[i],2],1]
      y2 = node_positions[connectivity[e_selected[i],2],2]
  
  
      if ratio_σ[e_selected[i]] < 0
        b = 1
        r = -ratio_σ[e_selected[i]]/ min_ratio_σ + 1
        g = -ratio_σ[e_selected[i]]/ min_ratio_σ + 1
      elseif ratio_σ[e_selected[i]] >= 0
        r = 1
        g = -ratio_σ[e_selected[i]]/ max_ratio_σ + 1
        b = -ratio_σ[e_selected[i]]/ max_ratio_σ + 1
      end
  
  
      ax2.plot([x1,x2], [y1,y2], linestyle = "-", color = (r,g,b), linewidth = 3);
    end
    
  
    xmax = maximum(node_positions[:,1]) -1
    ymax = maximum(node_positions[:,2]) -1
  
  
    line_z = collect(min_ratio_σ:0.001:max_ratio_σ) 
    stepy = ymax/(length(line_z))
    line_y = collect(0:stepy:ymax)
  
    for i = 1:1:size(line_z,1)-1
  
      if line_z[i] < 0
        b = 1
        r = -line_z[i]/ min_ratio_σ + 1
        g = -line_z[i]/ min_ratio_σ + 1
      elseif line_z[i] > 0
        r = 1
        g = -line_z[i]/ max_ratio_σ + 1
        b = -line_z[i]/ max_ratio_σ + 1
      end
  
    
      ax2.plot([xmax, xmax], [line_y[i], line_y[i+1]], linestyle = "-", color = (r,g,b), linewidth = 16);
    end
  
    ax2.annotate(round(min_ratio_σ, digits = 3),
              xy=[1;0],
              xycoords="axes fraction",
              xytext=[-30,5],
              textcoords="offset points",
              fontsize=10.0,
              ha="right",
              va="bottom",
              color = "white") 
  
    ax2.annotate(round(max_ratio_σ, digits = 3),
            xy=[1;1],
            xycoords="axes fraction",
            xytext=[-30,-5],
            textcoords="offset points",
            fontsize=10.0,
            ha="right",
            va="top",
            color = "white") 

    ax2.set_aspect("equal")


  end

  return σ

end




function plot_strain(node_positions, new_positions, connectivity)

  ϵ = zeros(size(connectivity,1))

   for i = 1:1:size(connectivity,1)

       x1 = node_positions[connectivity[i,1],1]
       y1 = node_positions[connectivity[i,1],2]

       x2 = node_positions[connectivity[i,2],1]
       y2 = node_positions[connectivity[i,2],2]

       L0 = sqrt((x2-x1)^2+(y2-y1)^2);

       x1 = new_positions[connectivity[i,1],1]
       y1 = new_positions[connectivity[i,1],2]

       x2 = new_positions[connectivity[i,2],1]
       y2 = new_positions[connectivity[i,2],2]

       L1 = sqrt((x2-x1)^2+(y2-y1)^2);

       ϵ[i] = (L1-L0) /L0;

   end

  
  (fig1, ax1) = subplots(1, 1, figsize=(12,12))
  #(fig1, ax1) = subplots(1, 1, figsize=(7,5)) # Figures papers
  ax1.set_facecolor((0, 0, 0))

  max_ϵ = maximum(ϵ);
  min_ϵ = minimum(ϵ);


  #max_ϵ = 1.1;
  #min_ϵ = -0.15;


  for i = 1:1:size(connectivity,1)

    x1 = new_positions[connectivity[i,1],1]
    y1 = new_positions[connectivity[i,1],2]

    x2 = new_positions[connectivity[i,2],1]
    y2 = new_positions[connectivity[i,2],2]

    if ϵ[i] >= max_ϵ
      ϵ[i] = max_ϵ
    elseif ϵ[i] <= min_ϵ
      ϵ[i] = min_ϵ
    end

    if ϵ[i] <= 0
      b = 1
      r = -ϵ[i]/ min_ϵ + 1
      g = -ϵ[i]/ min_ϵ + 1
    elseif ϵ[i] > 0
      r = 1
      g = -ϵ[i]/ max_ϵ + 1
      b = -ϵ[i]/ max_ϵ + 1
    end

    ax1.plot([x1,x2], [y1,y2], linestyle = "-", color = (r,g,b), linewidth = 3);
  end
  
  xmax = maximum(new_positions[:,1]) + 0.2
  ymax = maximum(new_positions[:,2]) + 0.2
  xmin = minimum(new_positions[:,1]) - 0.2
  ymin = minimum(new_positions[:,2]) - 0.2

  print("max x = ", xmax, "\n")
  print("max y = ", ymax, "\n")
  # xmax = 4.2;
  # ymax = 4.2;

  # xmin = -0.2;
  # ymin = -0.2;

  ax1.set_xlim([xmin, xmax])
  ax1.set_ylim([ymin, ymax])

  ax1.set_xticklabels([])
  ax1.set_xticks([])
  ax1.set_yticklabels([])
  ax1.set_yticks([])

  line_z = collect(min_ϵ:0.001:max_ϵ) 
  stepy = ymax/(length(line_z))
  line_y = collect(0:stepy:ymax)

  for i = 1:1:size(line_z,1)-1

    if line_z[i] <= 0
      b = 1
      r = -line_z[i]/ min_ϵ + 1
      g = -line_z[i]/ min_ϵ + 1
    elseif line_z[i] > 0
      r = 1
      g = -line_z[i]/ max_ϵ + 1
      b = -line_z[i]/ max_ϵ + 1
    end

  
    ax1.plot([xmax, xmax], [line_y[i], line_y[i+1]], linestyle = "-", color = (r,g,b), linewidth = 16);
  end

  annotate(round(min_ϵ, digits = 3),
            xy=[1;0],
            xycoords="axes fraction",
            xytext=[-30,5],
            textcoords="offset points",
            fontsize=10.0,
            ha="right",
            va="bottom",
            color = "white") 

  annotate(round(max_ϵ, digits = 3),
          xy=[1;1],
          xycoords="axes fraction",
          xytext=[-30,-5],
          textcoords="offset points",
          fontsize=10.0,
          ha="right",
          va="top",
          color = "white") 
  
  ax1.set_aspect("equal")
 
  return ϵ

end

