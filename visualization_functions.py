
# PREAMBLE:

import sys
import MDAnalysis
import numpy as np
from scipy import interpolate

def create_vis_state(visualization_frame_pdb,selection_list,paths,vis_state_file_name, node_sphere_radius = 1.0, node_sphere_rgb = (1.0,1.0,1.0), shortest_path_rgb = (1.0,1.0,1.0), longest_path_rgb = (1.0,1.0,1.0), shortest_path_radius = 1.00, longest_path_radius = 0.01, node_sphere_color_index = 34, VMD_color_index_range = (35,45), VMD_resolution = 200, VMD_spline_smoothness = 1.0):
        """
        unused colors in vmd: https://www.ks.uiuc.edu/Research/vmd/current/ug/node83.html#ug:topic:coloring
        """
        # ----------------------------------------
        # ANALYZE PATHS FOR PATH LENGTHS
        # ----------------------------------------
        if len(paths) != 1:
                # analyze paths
                lengths = [path[0] for path in paths]
                min_length = min(lengths)
                max_length = max(lengths)

                length_ratios = [(length - min_length)/(max_length - min_length) for length in lengths]
        else:
                length_ratios = [0.5]
        
        # ----------------------------------------
        # ANALYZE PATHS FOR NODE INDICES
        # ----------------------------------------
        nodes_used = []
        for path in paths:
                for index in path[1:]:
                        if not index in nodes_used:
                                nodes_used.append(index)
        
        node_resids = []
        node_positions = []
        for node in nodes_used:
                if str(type(selection_list[node])) == "<class 'MDAnalysis.core.groups.Atom'>":
                        node_resids.append(selection_list[node].resid)
                        node_positions.append(selection_list[node].position)
                else:
                        for resid in np.unique(selection_list[node].resids):
                                node_resids.append(int(resid))
                        node_positions.append(selection_list[node].center_of_mass())

        node_resids = list(set(node_resids))
        node_vmd_selection = 'resid'
        for node in node_resids:
                node_vmd_selection += ' ' + str(node)

        # ----------------------------------------
        # PREP A DICTIONARY FILLED WITH COLOR IDS AND RGB VALUES
        # ----------------------------------------
        color_defs = {}
        color_defs[node_sphere_color_index] = "color change rgb " + str(node_sphere_color_index) + " " + str(node_sphere_rgb[0]) + " " + str(node_sphere_rgb[1]) + " " + str(node_sphere_rgb[2]) + '\n'

        shortest_path_rgb = np.array(shortest_path_rgb)
        longest_path_rgb = np.array(longest_path_rgb)
        possible_colorids = range(VMD_color_index_range[0], VMD_color_index_range[1]+1)
        nColorid = len(possible_colorids)
        for colorid in possible_colorids:
                thiscolor = shortest_path_rgb + (colorid - VMD_color_index_range[0]) * (longest_path_rgb - shortest_path_rgb)/(nColorid-1)
                color_defs[colorid] = "color change rgb " + str(colorid) + " " + str(thiscolor[0]) + " " + str(thiscolor[1]) + " " + str(thiscolor[2]) + '\n'

        # ----------------------------------------
        # BEGIN WRITING THE VIS STATE FILE
        # ----------------------------------------
        with open(vis_state_file_name,'w') as W:
                W.write('mol new ' + visualization_frame_pdb + ' type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n')
                W.write('mol delrep 0 top\n')
                W.write('mol representation NewCartoon 0.160000 50.000000 4.100000 0\n')
                W.write('mol color Name\n')
                W.write('mol selection {all}\n')
                W.write('mol material Opaque\n')
                W.write('mol addrep top\n')
                W.write('mol representation Licorice 0.150000 85.000000 85.000000\n')
                W.write('mol color Name\n')
                W.write('mol selection {' + node_vmd_selection + '}\n')
                W.write('mol material Opaque\n')
                W.write('mol addrep top\n')
                
                # ----------------------------------------
                # DRAW NODE SPHERES
                # ----------------------------------------
                if node_sphere_radius != 0.0:
                        # draw spheres at node positions
                        W.write('\n# Draw spheres at the nodes\n')
                        W.write(color_defs[node_sphere_color_index])
                        W.write('graphics top color ' + str(node_sphere_color_index) + '\n')
                        for node in node_positions:
                                W.write('draw sphere {' + str(node[0]) + ' ' + str(node[1]) + ' ' + str(node[2]) + '} resolution ' + str(VMD_resolution) + ' radius ' + str(node_sphere_radius) + '\n')

                W.write('set wisp_num_paths ' + str(len(paths)) + '\n')

                # ----------------------------------------
                # DRAW PATHWAYS USING INTERPOLATION
                # ----------------------------------------
                for i in range(len(paths)):
                        W.write('\n# Draw a new path\n# \tLength: ' + str(paths[i][0]) + '\n# \t Nodes: ' + str(paths[i][1:]) + '\n\n')
                        ratio = length_ratios[i]
                        color = int(np.floor(ratio * (nColorid-1)) + VMD_color_index_range[0])
                        radius = ratio * (longest_path_radius - shortest_path_radius) + shortest_path_radius

                        W.write(color_defs[color] + '\n')
                        W.write('graphics top color ' + str(color) + '\n')
                        
                        x_vals = []
                        y_vals = []
                        z_vals = []
                        for node in paths[i][1:]:
                                if str(type(selection_list[node])) == "<class 'MDAnalysis.core.groups.Atom'>":
                                        pos = selection_list[node].position
                                else:
                                        pos = selection_list[node].center_of_mass()
                                
                                x_vals.append(pos[0])
                                y_vals.append(pos[1])
                                z_vals.append(pos[2])
                        
                        try:
                                degree = len(x_vals)
                                if degree > 3: 
                                        degree = 3

                                tck, u = interpolate.splprep([x_vals,y_vals,z_vals],s=0,k=degree)
                                unew = np.arange(0,1.01,VMD_spline_smoothness)
                                out = interpolate.splev(unew,tck)

                                for t in range(len(out[0])-1):
                                        x1 = str(out[0][t])
                                        y1 = str(out[1][t])
                                        z1 = str(out[2][t])

                                        x2 = str(out[0][t+1])
                                        y2 = str(out[1][t+1])
                                        z2 = str(out[2][t+1])

                                        W.write('draw cylinder {' + x1 + ' ' + y1 + ' ' + z1 + '} {' + x2 + ' ' + y2 + ' ' + z2 + '} radius ' + str(radius) + ' resolution ' + str(VMD_resolution) + ' filled 0\n')
                        except:
                                W.write('draw cylinder {' + str(x_vals[0]) + ' ' + str(y_vals[0]) + ' ' + str(z_vals[0]) + '} {' + str(x_vals[-1]) + ' ' + str(y_vals[-1]) + ' ' + str(z_vals[-1]) + '} radius ' + str(radius) + ' resolution ' + str(VMD_resolution) + ' filled 0\n')

