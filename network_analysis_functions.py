
# PREAMBLE:

import sys
import copy
import numpy as np
from numpy.linalg import *
import networkx

######################### To Identify Paths ##############################

def get_paths(correlation_matrix,sources,sinks,number_of_paths):
        """
        """

        nNodes = len(correlation_matrix) 
        nNodes_range = range(nNodes)

        #populate graph nodes and weighted edges
        G = networkx.Graph(data = correlation_matrix)

        # ----------------------------------------
        # SHORTEST PATH CALCULATION
        # ----------------------------------------
        shortest_length = 99999999.999
        shortest_path = []
        for source in sources:
                for sink in sinks:
                        if source != sink:
                                short_path = networkx.dijkstra_path(G,source,sink,weight='weight')
                                length = 0
                                for t in range(len(short_path)-1):
                                        length += correlation_matrix[short_path[t],short_path[t+1]]
                                if length < shortest_length:
                                        shortest_length = length
                                        shortest_path = short_path
        
        print 'Shortest pathway information:', shortest_length, shortest_path
        
        # ----------------------------------------
        # SUBOPTIMAL PATHS CALCULATION
        # ----------------------------------------
        num_paths = 1
        path = [shortest_length]
        path.extend(shortest_path)
        paths = [path]

        cutoff = shortest_length

        cutoff_yields_max_num_paths_below_target = 0 
        cutoff_yields_min_num_paths_above_target = 1000000.0 

        # ----------------------------------------
        # first step, keep incrementing a little until you have more than the desired number of paths
        # ----------------------------------------
        while num_paths < number_of_paths:
                cutoff_in_array = np.array([cutoff])
                print 'new cutoff = ', cutoff
                
                # ----------------------------------------
                # starting get_paths_between_multiple_endpoints
                # ----------------------------------------
                temp_paths = []
                for source in sources:
                        for sink in sinks:
                                if source != sink:
                                        # ----------------------------------------
                                        # starting get_paths_fixed_endpoints
                                        # ----------------------------------------
                                        # measure shortest paths between source/sink nodes to all other nodes; based on the G graph
                                        source_lengths, source_paths = networkx.single_source_dijkstra(G,source,target=None,cutoff=None,weight='weight')
                                        sink_lengths, sink_paths = networkx.single_source_dijkstra(G,sink,target=None,cutoff=None,weight='weight')

                                        # converting networkx results from dictionaries to lists
                                        source_lengths_list = [source_lengths[key] for key in source_lengths.keys()]
                                        source_paths_list = [source_paths[key] for key in source_paths.keys()]
                                        sink_lengths_list = [sink_lengths[key] for key in sink_lengths.keys()]
                                        sink_paths_list = [sink_paths[key] for key in sink_paths.keys()]

                                        # collecting all node indices found by the shortest-path calculations done above; should return two lists filled with every node index; why is  this important to do?
                                        check_list_1 = []
                                        check_list_2 = []
                                        for i in nNodes_range:
                                                check_list_1.extend([source_paths_list[i][-1]])
                                                check_list_2.extend([sink_paths_list[i][-1]])

                                        # identifying nodes that should be included in the suboptimal pathways analysis
                                        node_list = []
                                        if check_list_1 == check_list_2:
                                                for i in nNodes_range:      # node i is the forced node; see if the sum of lengths between source and node i AND between sink and node i is less than the cutoff; if so, add the two paths (source to node i AND sink to node i) to the node list
                                                        if source_lengths_list[i] + sink_lengths_list[i] <= cutoff:
                                                                node_list.extend(source_paths_list[i][:])
                                                                node_list.extend(sink_paths_list[i][:])
                                        else:
                                                print 'Paths do not match up. Something is wrong here...'
                                                sys.exit()

                                        unique_nodes = list(set(node_list))
                                        unique_nodes.sort()
                                        
                                        # simplifying our original correlation matrix to only include correlations between nodes that were identified above
                                        node_length_range = range(len(unique_nodes))
                                        new_matrix = np.zeros((nNodes,nNodes))
                                        for i in node_length_range:
                                                for j in node_length_range:
                                                        new_matrix[unique_nodes[i]][unique_nodes[j]] = correlation_matrix[unique_nodes[i]][unique_nodes[j]]

                                        temp_G = networkx.Graph(data=new_matrix,labels=unique_nodes)

                                        # begin to study the suboptimal paths by building the pathways and calculating their lengths.
                                        length = 0.0
                                        paths_growing_out_from_source = [[length,source]]
                                        full_paths_from_source_to_sink = []
                                        
                                        while len(paths_growing_out_from_source) > 0:
                                                # ----------------------------------------
                                                # starting expand_growing_paths_one_step
                                                # ----------------------------------------
                                                for i in range(len(paths_growing_out_from_source)):
                                                        if paths_growing_out_from_source[i][0] > cutoff:
                                                                paths_growing_out_from_source.pop(i)    # remove this path from the list
                                                                break
                                                        elif paths_growing_out_from_source[i][-1] == sink:
                                                                full_paths_from_source_to_sink.append(paths_growing_out_from_source[i])
                                                                paths_growing_out_from_source.pop(i)    # remove this path from the list
                                                                break
                                                        elif paths_growing_out_from_source[i][-1] != sink:
                                                                node_neighbors = temp_G.neighbors(paths_growing_out_from_source[i][-1])
                                                                for j in range(len(node_neighbors)):
                                                                        if not node_neighbors[j] in paths_growing_out_from_source[i]:
                                                                                temp = paths_growing_out_from_source[i][:]
                                                                                temp.append(node_neighbors[j])
                                                                                temp[0] += temp_G.edge[temp[-2]][temp[-1]]['weight']
                                                                                paths_growing_out_from_source.insert((i+j+1),temp)
                                                                paths_growing_out_from_source.pop(i)    # remove this path from the list
                                                                break
                                                        else:
                                                                print 'Something is wrong.'
                                                                sys.exit()

                                                # ----------------------------------------
                                                # ending expand_growing_paths_one_step
                                                # ----------------------------------------
                                        
                                        full_paths_from_source_to_sink.sort()
                                        temp_paths.extend([full_paths_from_source_to_sink[i] for i in range(len(full_paths_from_source_to_sink))])
                                        print 'Done with path calc for this cutoff'
                                        # ----------------------------------------
                                        # ending get_paths_fixed_endpoints
                                        # ----------------------------------------
                # ----------------------------------------
                # ending get_paths_between_multiple_endpoints
                # ----------------------------------------
                # ----------------------------------------
                # starting remove_redundant_paths
                # ----------------------------------------
                if len(temp_paths) != 1:
                        for i in range(len(temp_paths)-1):
                                path1 = temp_paths[i]
                                if not path1 is None:
                                        for j in range(i+1,len(temp_paths)):
                                                path2 = temp_paths[j]
                                                if not path2 is None:
                                                        if len(path1) == len(path2):
                                                                pth1 = copy.deepcopy(path1[1:])
                                                                pth2 = copy.deepcopy(path2[1:])
                                                                if pth1[0] < path1[-1]:
                                                                        pth1.reverse()
                                                                if pth2[0] < path2[-1]:
                                                                        pth2.reverse()
                                                                if pth1 == pth2:
                                                                        temp_paths[j] = None
                        while None in temp_paths:
                                temp_paths.remove(None)

                # ----------------------------------------
                # ending remove_redundant_paths
                # ----------------------------------------
                num_paths = len(temp_paths)
                paths.extend(temp_paths)
                
                if num_paths < number_of_paths and cutoff > cutoff_yields_max_num_paths_below_target:
                        cutoff_yields_max_num_paths_below_target = cutoff 
                if num_paths > number_of_paths and cutoff < cutoff_yields_min_num_paths_above_target:
                        cutoff_yields_min_num_paths_above_target = cutoff

                cutoff += 0.1

        if len(paths) != 1:
                for i in range(len(paths)-1):
                        path1 = paths[i]
                        if not path1 is None:
                                for j in range(i+1,len(paths)):
                                        path2 = paths[j]
                                        if not path2 is None:
                                                if len(path1) == len(path2):
                                                        pth1 = copy.deepcopy(path1[1:])
                                                        pth2 = copy.deepcopy(path2[1:])
                                                        if pth1[0] < path1[-1]:
                                                                pth1.reverse()
                                                        if pth2[0] < path2[-1]:
                                                                pth2.reverse()
                                                        if pth1 == pth2:
                                                                paths[j] = None
                                    
        while None in paths:
                paths.remove(None)

        paths.sort() # sort the paths by length

        if num_paths != number_of_paths: # so further refinement is needed
                paths = paths[:number_of_paths]

        return paths

