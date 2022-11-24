"""
Model development: Dr Alessandra Bonfanti, Prof Alexandre Kabla
Code author: Dr Alessandra Bonfanti

To gain a conceptual understanding of how tension in mitotic cells evolves relative to the tension 
in interphase cells and in the tissue, we devised a simple computational model of the monolayer as 
a purely elastic 2D material in which cells have a preset active tension and rigidity. By varying tension, 
rigidity and position of the boundaries, we aim to reproduce a range of experimental conditions and 
characterize the distribution of stress in the viscinity of the mitotic cell. 
The simulated scenarios are:
    1) monolayer in its initial configuration clamped at both ends is subjected to a tensile stress due to cell contractility
    2) the monoalyer with internal cell contractility is compressed
    3) the internal contractility is increased while the monolayer is under a compressive state
    4) the internal contractility of a clamped monolayer is reduced
    5) the monolayer with reduced contractility is stretched

Scientific Publication:
"Tension at intercellular junctions is necessary for accurate orientation of cell division in the epithelium plane" 
"""


using PyPlot
using LinearAlgebra
using Colors
using ColorSchemes
using DelimitedFiles
include("functions.jl")



function main_triangle_inclusion_fixed(param1,param2, param3, filename)


    EA = param1
    disp = param2
    prestress = param3;
    EAmult = 1              # Stiffness properties of the dividing cell (if = 1, then it has the same properties of the serrounding cells)
    preMult = 1             # Prestress of the dividing cell (if = 1, then it has the same prestress of the serrounding cells)
    save_stress = false
    save_tensor = false
    save_strain = false
    name_file_stress = string(filename, "_stress.txt")    # Saving stress in 
    name_file_strain = string(filename, "_strain.txt")
    name_file_tensor = string(filename, "_tensor.txt")

    # Generate geometry
    nodes_position, nodes_q, connectivity = PositionParticles(41, 41, 0.1, 0.1, reg = true, pert = 0.6)
    print("@ Lattice generated", "\n")

    # # Select nodes and elements of the inclusion (= dividing cell)
    selectedcell, e_selected = SelectCell(nodes_position, connectivity, 1.9, 2.1, 1.9, 2.1)

    # Plot lattice geometry
    plot_lattice(nodes_position, connectivity, "white", "--", selected = e_selected)
    # Plot the dividing cell in the monolayer with a different color
    plot(nodes_position[selectedcell,1], nodes_position[selectedcell,2], "o", color = "green")


    # Create global stiffness matrix
    K, F_vect = K_global(nodes_position, nodes_q, connectivity, EA, prestress; selected = e_selected, EAmultiplier = EAmult, preMultiplier = preMult)
    print("@ Global stiffness matrix generated", "\n")

    # Select nodes to which apply the force and create force vector
    x_max = maximum(nodes_position[:,1])
    y_max = maximum(nodes_position[:,2])

    # Apply BC
    node_left = SelectNodes(nodes_position, x = 0)
    node_right = SelectNodes(nodes_position, x = x_max)
    node_middle = SelectNodes(nodes_position, x = x_max/2)
    node_middle_length = SelectNodes(nodes_position, y = y_max/2)


    K_BC, F_BC, q_to_update, q_disp_right, q_disp_left = applyBC(K, F_vect, node_right, node_left, nodes_q, disp)
    print("@ BC applied", "\n")

    # Solve the systme
    q = K_BC \ F_BC;
    print("@ System solved", "\n")

    # Plot new positions
    new_positions, all_q = update_position(nodes_position, nodes_q, q, q_to_update, q_disp_right, q_disp_left, disp); 
    print("@ Updated geometry", "\n")

    stress = plot_stress(nodes_position, new_positions, EA, connectivity, prestress; selected = e_selected, EAmultiplier = EAmult, preMultiplier = preMult)
    plot(new_positions[selectedcell,1], new_positions[selectedcell,2], "s", color = "black", markersize = 3)


    strain = plot_strain(nodes_position, new_positions, connectivity)
    plot(new_positions[selectedcell,1], new_positions[selectedcell,2], "s", color = "black", markersize = 3)


    if save_stress == true
        writedlm( name_file_stress,  stress)
    end

    if save_strain == true
        writedlm( name_file_strain,  strain)
    end


    tensors = Array{Float64}(undef, 0, 4);

    #plot_lattice(nodes_position, connectivity, "white", "--", selected = e_selected)

    print("@ Stress tensor @ dividing cell", "\n")
    # Select elements for tensor calculation
    selectedcell, e_selected = SelectCell(nodes_position, connectivity, 1.99, 2.01, 1.99, 2.01)
    #plot(nodes_position[selectedcell,1], nodes_position[selectedcell,2], "o", color = "green")
    tensor = calculate_tensor(selectedcell, new_positions, connectivity, stress)
    calculate_strain(selectedcell, nodes_position, connectivity, strain)

    tensors = vcat(tensors, [tensor[1,1] tensor[1,2] tensor[2,1] tensor[2,2]])

    print("@ Stress tensor @ LEFT dividing cell", "\n")
    # Select elements for tensor calculation
    selectedcell, e_selected = SelectCell(nodes_position, connectivity, 1.79, 1.81, 1.99, 2.01)
    #plot(nodes_position[selectedcell,1], nodes_position[selectedcell,2], "o", color = "green")
    tensor = calculate_tensor(selectedcell, new_positions, connectivity, stress)
    tensors = vcat(tensors, [tensor[1,1] tensor[1,2] tensor[2,1] tensor[2,2]])

    print("@ Stress tensor @ RIGHT dividing cell", "\n")
    # Select elements for tensor calculation
    selectedcell, e_selected = SelectCell(nodes_position, connectivity, 2.19, 2.21, 1.99, 2.01)
    #plot(nodes_position[selectedcell,1], nodes_position[selectedcell,2], "o", color = "green")
    tensor = calculate_tensor(selectedcell, new_positions, connectivity, stress)
    tensors = vcat(tensors, [tensor[1,1] tensor[1,2] tensor[2,1] tensor[2,2]])

    print("@ Stress tensor @ TOP dividing cell", "\n")
    # Select elements for tensor calculation
    selectedcell, e_selected = SelectCell(nodes_position, connectivity, 1.99, 2.01, 2.19, 2.21)
    #plot(nodes_position[selectedcell,1], nodes_position[selectedcell,2], "o", color = "green")
    tensor = calculate_tensor(selectedcell, new_positions, connectivity, stress)
    tensors = vcat(tensors, [tensor[1,1] tensor[1,2] tensor[2,1] tensor[2,2]])

    print("@ Stress tensor @ BOTTOM dividing cell", "\n")
    # Select elements for tensor calculation
    selectedcell, e_selected = SelectCell(nodes_position, connectivity, 1.99, 2.01, 1.79, 1.81)
    #plot(nodes_position[selectedcell,1], nodes_position[selectedcell,2], "o", color = "green")
    tensor = calculate_tensor(selectedcell, new_positions, connectivity, stress)
    tensors = vcat(tensors, [tensor[1,1] tensor[1,2] tensor[2,1] tensor[2,2]])

    print("@ Stress tensor @ LEFT INFTY dividing cell", "\n")
    # Select elements for tensor calculation
    selectedcell, e_selected = SelectCell(nodes_position, connectivity, 0.99, 1.01, 1.99, 2.01)
    #plot(nodes_position[selectedcell,1], nodes_position[selectedcell,2], "o", color = "green")
    tensor = calculate_tensor(selectedcell, new_positions, connectivity, stress)
    tensors = vcat(tensors, [tensor[1,1] tensor[1,2] tensor[2,1] tensor[2,2]])

    print("@ Stress tensor @ RIGHT INFTY dividing cell", "\n")
    # Select elements for tensor calculation
    selectedcell, e_selected = SelectCell(nodes_position, connectivity, 2.99, 3.01, 1.99, 2.01)
    #plot(nodes_position[selectedcell,1], nodes_position[selectedcell,2], "o", color = "green")
    tensor = calculate_tensor(selectedcell, new_positions, connectivity, stress)
    tensors = vcat(tensors, [tensor[1,1] tensor[1,2] tensor[2,1] tensor[2,2]])

    print("@ Stress tensor @ TOP INFTY dividing cell", "\n")
    # Select elements for tensor calculation
    selectedcell, e_selected = SelectCell(nodes_position, connectivity, 1.99, 2.01, 2.99, 3.01)
    #plot(nodes_position[selectedcell,1], nodes_position[selectedcell,2], "o", color = "green")
    tensor = calculate_tensor(selectedcell, new_positions, connectivity, stress)
    tensors = vcat(tensors, [tensor[1,1] tensor[1,2] tensor[2,1] tensor[2,2]])

    print("@ Stress tensor @ BOTTOM INFTY dividing cell", "\n")
    # Select elements for tensor calculation
    selectedcell, e_selected = SelectCell(nodes_position, connectivity, 1.99, 2.01, 0.99, 1.01)
    #plot(nodes_position[selectedcell,1], nodes_position[selectedcell,2], "o", color = "green")
    tensor = calculate_tensor(selectedcell, new_positions, connectivity, stress);
    tensors = vcat(tensors, [tensor[1,1] tensor[1,2] tensor[2,1] tensor[2,2]])

    if save_tensor == true
        writedlm( name_file_tensor, tensors)
    end


end




# ---------------------------------- Main ------------------------------------------------------------------

n = 44    # Number of steps for the simulation

for i = 1:1:n
    close("all")

    # Case 01: No axial displacement, only prestress
     param1 =  1                    # EA = where E is the Young's modulus and A is the corss section area
     param2 = 0.0                   # Displacement at the two edges (Boundary conditions)
     param3 = (i-1) * 0.35/(n-1)    # Prestress

    # Case 02: constant prestress, compression
    #  param1 =  1    
    #  param2 = -(i-1) * 0.6/(n-1)
    #  param3 = 0.35

    # Case 03: increasing prestress, compression
    #  param1 =  1 + (i-1) *0.5/(n-1);   
    #  param2 = -0.6
    #  param3 = 0.35 + (i-1) * 0.25/(n-1)

    # Case 04: 
    #  param1 =  1 - (i-1) *0.3/(n-1);    
    #  param2 = 0.0
    #  param3 = 0.35 - (i-1) * 0.25/(n-1)

    #  param1 =  0.7   # case 02
    #  param2 = (i-1) * 1/(n-1)
    #  param3 = 0.1


    filename = string("./case01_Inclusion/case01Prova_",i)   # Name of the file where to save the stress tensor at the 9 selected positions 
                                                                 #  @ dividing cell
                                                                 #  @ top-bottom-left-right next to the dividing cell
                                                                 #  @ top-bottom-left-right far from the dividing cell
    
    # Run the simulation
    main_triangle_inclusion_fixed(param1, param2, param3, filename) 

    # Save the strain plot in pdf file
    figname = string("./case01_Inclusion/FigureplotStrainProva_", i, ".png")  
    savefig(figname)

end
