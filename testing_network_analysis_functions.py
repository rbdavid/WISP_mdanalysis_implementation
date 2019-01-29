
# PREAMBLE:

import sys
import copy
import numpy as np
from numpy.linalg import *
import networkx
import time

######################### To Identify Paths ##############################

def get_paths(correlation_matrix,sources,sinks,number_of_paths,user_defined_percent_increase_cutoff=0.001):
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
                # killing job if instances a node is found in both the source and sink lists.
                if source in sinks:    
                        print 'A source node (%d) is also in the sink node list. This will return trivial results and is not interesting. Killing Job...' %(source)
                        sys.exit()

                # measuring shortest pathway length and nodes
                for sink in sinks:
                        short_path = networkx.dijkstra_path(G,source,sink,weight='weight')
                        length = 0
                        for t in range(len(short_path)-1):
                                length += correlation_matrix[short_path[t],short_path[t+1]]
                        if length < shortest_length:
                                shortest_length = length
                                shortest_path = short_path
        
        print 'Shortest pathway information:', shortest_length, shortest_path, '\n'
      
        # ----------------------------------------
        # INITIALIZE THE PATHS LIST, FILLED WITH THE OPTIMAL PATH
        # ----------------------------------------
        num_paths = 1
        path = [shortest_length]    # prepping the first sublist to be placed in the paths list of lists; appending shortest pathway length to the path list
        path.extend(shortest_path)  # extending the path list to include the node indices for the respective path
        paths = [path]      # paths is a list of lists, where each sublist holds the length and node indices for the respective paths

        # ----------------------------------------
        # LONGEST PATH CALCULATION  ### the current algorithm used in this section is bs... will likely always find the direct source-sink path as the longest path.
        # ----------------------------------------
        temp_G = networkx.Graph(data=(-correlation_matrix+np.max(correlation_matrix)))  # edge weights are the negative of their original values; :. edges with the long distances are short; use dijkstra algorithm to find the shortest path in this upside-down representation of our original graph 
        longest_length = 0.0
        longest_path = []
        for source in sources:
                for sink in sinks:
                        long_path = networkx.dijkstra_path(temp_G,source,sink,weight='weight')
                        temp_length = 0
                        for t in range(len(long_path)-1):
                                temp_length += correlation_matrix[long_path[t],long_path[t+1]]  # sum up the right-side up correlation_matrix values to calculate the true distance for the longest path
                        if temp_length > longest_length:
                                longest_length = temp_length
                                longest_path = long_path

        print 'Longest pathway information:', longest_length, longest_path

        # ----------------------------------------
        # SUBOPTIMAL PATHS CALCULATION
        # ----------------------------------------
        cutoff = shortest_length
        cutoff_iteration = (longest_length - shortest_length)*user_defined_percent_increase_cutoff

        print cutoff_iteration
       
        temp_file = open('temp.dat','w',1)

        # ----------------------------------------
        # prepping node lists to search for pathways through... this was originally being done within each iteration of the while loop but should only need to be done once!
        # ----------------------------------------
        
        # measure shortest paths between source/sink nodes to all other nodes; based on the original G graph
        source_lengths, source_paths = networkx.single_source_dijkstra(G,source,target=None,cutoff=None,weight='weight')
        sink_lengths, sink_paths = networkx.single_source_dijkstra(G,sink,target=None,cutoff=None,weight='weight')

        # converting networkx results from dictionaries to lists
        source_lengths_list = [source_lengths[key] for key in source_lengths.keys()]
        source_paths_list = [source_paths[key] for key in source_paths.keys()]
        sink_lengths_list = [sink_lengths[key] for key in sink_lengths.keys()]
        sink_paths_list = [sink_paths[key] for key in sink_paths.keys()]

        # ----------------------------------------
        # first step, keep incrementing a little until you have more than the desired number of paths
        # ----------------------------------------
        while num_paths < number_of_paths:
                start = time.time()
                print 'new cutoff = ', cutoff
                
                # ----------------------------------------
                # starting get_paths_between_multiple_endpoints
                # ----------------------------------------
                temp_paths = []
                for source in sources:
                        for sink in sinks:
                                # identifying nodes that should be included in the suboptimal pathways analysis for the given cutoff
                                node_list = []
                                for i in nNodes_range:      # assumes all nodes are potential pathway nodes. node i is the forced node; see if the sum of lengths between source and node i AND between sink and node i is less than the cutoff; if so, add the nodes observed in both paths (source to node i AND sink to node i) to the node list
                                        if source_lengths_list[i] + sink_lengths_list[i] <= cutoff:
                                                node_list.extend(source_paths_list[i][:])
                                                node_list.extend(sink_paths_list[i][:])

                                unique_nodes = list(set(node_list))
                                unique_nodes.sort()
                                
                                # simplifying our original correlation matrix to only include correlations between nodes that were identified above
                                new_matrix = np.zeros((nNodes,nNodes))
                                for i in unique_nodes:
                                        for j in unique_nodes:
                                                new_matrix[i][j] = correlation_matrix[i][j]

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
                                                if paths_growing_out_from_source[i][0] > cutoff:    # if current path has length greater than cutoff, remove this pathway from the paths_growing_out_from_source list
                                                        paths_growing_out_from_source.pop(i)
                                                        break

                                                elif paths_growing_out_from_source[i][-1] == sink:  # if the last (aka latest) node in the path is the sink node, add this path to the full_paths_from_source_to_sink list and remove from the paths_growing_out_from_source list
                                                        full_paths_from_source_to_sink.append(paths_growing_out_from_source[i])
                                                        paths_growing_out_from_source.pop(i)
                                                        break

                                                elif paths_growing_out_from_source[i][-1] != sink:  # if the last (aka latest) node in the path is not the sink node, continue adding nodes to this pathway, then remove the old/originating path
                                                        node_neighbors = temp_G.neighbors(paths_growing_out_from_source[i][-1])     # find "neighbors" to this latest node... ### look into this function; is this just returning every single node? how does the networkx algorithm determine a neighbor?
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
                                # ----------------------------------------
                                # ending get_paths_fixed_endpoints
                                # ----------------------------------------
                print 'Done with path calc for this cutoff'
                end = time.time()
                pathway_time = end - start
                # ----------------------------------------
                # ending get_paths_between_multiple_endpoints
                # ----------------------------------------
                # ----------------------------------------
                # starting remove_redundant_paths
                # ----------------------------------------
                start = time.time()
                
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

                end = time.time()
                redundant_paths_time = end - start

                cutoff += cutoff_iteration * (1-float(num_paths)/number_of_paths)
                temp_file.write('%f   %d   %f   %f   %f\n' %(cutoff,num_paths,1-float(num_paths)/number_of_paths,pathway_time,redundant_paths_time))

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

        temp_file.close()

        return paths

