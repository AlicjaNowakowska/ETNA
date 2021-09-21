import warnings
from IPython.utils import io
warnings.filterwarnings("ignore")
from graph_tool.all import *
import graph_tool.all as gt
import ipywidgets as widgets
import matplotlib.pyplot as plt
import numpy as np
from IPython.display import clear_output
#from ipywidgets import Layout
import seaborn as sns
#from ipywidgets import HBox, VBox
from ipywidgets import *
import rpy2.robjects.packages as rpackages
from rpy2.robjects.packages import importr
import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.vectors import StrVector

utils = rpackages.importr('utils')
with io.capture_output() as captured:
  utils.install_packages('poweRlaw', repos="https://cloud.r-project.org")
x = rpackages.importr('poweRlaw')

class My_Network:

  def __init__(self, file_name):
    # Graph is created through the file upload
    self.G = graph_tool.load_graph(file_name)

  def prepare_the_network(self):
    """

    Graph preparation includes:
    1) Making it undirected
    2) Removal of parallel edges if they are present (if not nothing happens) 
    3) Extraction of the largest connected component that is treated as the final ready-to-use graph (all other components are removed).
    
    """

    self.G.set_directed(False) # 1)
    graph_tool.stats.remove_parallel_edges(self.G) # 2)

    #3):
    comp, hist = graph_tool.topology.label_components(self.G) 
    label = gt.label_largest_component(self.G)
    to_remove = []
    for v in self.G.vertices():
      if label[v]==0:
        to_remove.append(v)
    for v in reversed(sorted(to_remove)):
      self.G.remove_vertex(v)

  """ 
  
  The following functions are responsible for calculation of centrality measures and clustering coefficient. 
  It is done by generating a corresponding map of the form: node - value of the measure. 

  """

  def create_degree_distribution_map(self):
    my_map = self.G.degree_property_map("total")
    return my_map

  def create_betweenness_distribution_map(self):
    v_betweeness_map, e_betweenness_map = graph_tool.centrality.betweenness(self.G)
    my_map = v_betweeness_map
    return my_map

  def create_closeness_distribution_map(self):
    my_map = graph_tool.centrality.closeness(self.G)
    return my_map 
  
  def create_eigenvector_distribution_map(self):
    eigen_value, v_eigen_map = graph_tool.centrality.eigenvector(self.G)
    my_map = v_eigen_map
    return my_map

  def create_clustering_map(self):
    my_map = graph_tool.clustering.local_clustering(self.G)
    return my_map

  def plot_map_histogram(self, my_map, measure_name, block = True):
    """
    
    Plotting map histogram function contains a code for the graph generation
    using matplotlib library given the graph-tool map for the measure of interest.

    """

    # General settings:
    plt.style.use('seaborn-whitegrid')  
    fig, ax = plt.subplots(constrained_layout=True, figsize=(5, 5))
    FONT = 15

    # Preparing the data:
    my_map = my_map.fa # Extraction of the map's values - now the normal pythonic list is obtained as the representation of the measure's values.

    # Calculating basic statistics:
    to_calculate_statistics = list(my_map)
    avg = round(np.mean(to_calculate_statistics),4)
    std = round(np.std(to_calculate_statistics),2)

    # Creating the histogram:
    n=15
    a = ax.hist(my_map, bins=n, facecolor="lightblue",weights=np.zeros_like(my_map) + 1. / len(my_map))
    bins_mean = [0.5 * (a[1][j] + a[1][j+1]) for j in range(n)]
    sticks_to_mark = ([], [])
    for k in range(len(a[0])):
      if a[0][k] == 0:
        pass
      else:
        sticks_to_mark[0].append(bins_mean[k])
        sticks_to_mark[1].append(a[0][k])
    ax.plot(sticks_to_mark[0], sticks_to_mark[1], "b+") 
    ax.set_xlabel("Value", fontsize = FONT)
    ax.set_ylabel("Fraction of nodes", fontsize = FONT)
    ax.set_title(measure_name +" histogram \n Mean value: " + str(avg)+ ", Std: "+ str(std), fontsize = FONT)
    plt.show(block=block)
    return fig, ax

  def hubs_impact_check(self):
    """ 
    
    hubs_impact_check function is used for the evaluation of hubs and low-degree nodes' contribution to the number of links present in the graph.
    This is done by extracting all the possible values of the degree (1) and then looping over them (2). Within the loop for each degree number
    all nodes with the degree below or equal to it are extracted to form the subnetwork (3). The number of links and nodes in the subnetwork
    is divided by the corresponding total numbers in the network (4) to evaluate the contribution of the following degree groups. 

    """

    largest_N = self.G.num_vertices()
    largest_E = self.G.num_edges()
    degrees = self.G.get_total_degrees(self.G.get_vertices()) 
    Ns = []
    Es = []
    degrees_set = list(set(degrees)) # 1)
    degrees_set.sort() 
    degrees_map = self.G.degree_property_map("total")

    for degree in degrees_set: # 2)
      cut = degree  
      u = gt.GraphView(self.G, vfilt = lambda v: degrees_map[v]<=cut) # 3)
      current_N = u.num_vertices()/largest_N
      current_E = u.num_edges()/largest_E # 4)
      Ns.append(current_N)
      Es.append(current_E)

    for i in Ns: # to jest ten indeks gdzie jest przjście 90%
      if i >=0.9:
        index = Ns.index(i)
        break

    return Ns, Es, degrees_set
  
  def plot_hubs_impact1(self, degrees_set, Es, block = True): #to use it first need to execute hubs_impact_check
    """ 

    Plot_hubs_impact1 requires data that is generated by hubs_impact_check function. 
    It generates the plot that represents how the following degree groups contribute to the number of links present in the whole network.

    """

    # Plot settings:
    FONT = 15
    plt.style.use('seaborn-whitegrid')  
    plt.figure(figsize=(5,5))
    plt.xticks(fontsize=FONT-3)
    plt.yticks(fontsize=FONT-3)
    plt.xlabel("K", fontsize= FONT)
    plt.ylabel("Fraction of links", fontsize= FONT)
    plt.title("Size of the subnetwork of all nodes with\n a degree smaller or equal K", fontsize= FONT)
    # Plotting the data
    plt.plot(degrees_set, Es, "o", markersize=4, color="royalblue")
    plt.show(block = block)
    
  def plot_hubs_impact2(self, degrees_set, Es, Ns, block = True):
    """ 

    Plot_hubs_impact2 requires data that is generated by hubs_impact_check function. 
    It generates the plot that represents how the following percentages of the total number of nodes contribute to 
    the total number of links present in the whole network.

    """

    # Plot settings:
    FONT=15
    plt.style.use('seaborn-whitegrid')  
    plt.figure(figsize=(5,5))
    sns.set_context("paper", rc={"font.size":FONT,"axes.titlesize":FONT,"axes.labelsize":FONT, "xtick.labelsize":FONT-3, "ytick.labelsize":FONT-3,
                      "legend.fontsize":FONT-3, "legend.titlesize":FONT-3})  
    # Plotting the data
    fig = sns.scatterplot(x= Ns, y=Es, hue=np.log(degrees_set), palette="dark:blue_r")  
    fig.set(xlabel='Fraction of nodes', ylabel='Fraction of links', title="Size of the subnetwork of all nodes with\n a degree smaller or equal K")
    plt.legend(title="Log(K)", loc ="upper left", title_fontsize=FONT-3) 
    plt.show(block = block)
  
  def calculate_assortativity_value(self):
    # Calculation of the degree correlation coefficient:
    return gt.assortativity(self.G, "total")[0]
  
  def plot_ANND(self, normed = False, errorbar = True, block = True):
    """

    plot_ANND generates Average Nearest Neighbour Degree plot that represents the mixing patterns between different groups of the nodes. 
    Each group consists of the nodes of the same degree.

    """

    # Plot settings: 
    FONT = 15
    plt.style.use('seaborn-whitegrid')  
    fig = plt.figure(figsize=(5,5))
    plt.xlabel("Source degree (k)", fontsize = FONT)
    plt.ylabel("$<k_{nn}(k)>$", fontsize = FONT)
    title = "Average degree of\n the nearest neighbours" if normed == False else "Normed average degree of\n the nearest neighbours"
    plt.title(title, fontsize = FONT)

    # Calculating correlation vectors for ANND plot
    h = gt.avg_neighbor_corr(self.G, "total", "total")
    x = h[2][:-1]
    y = h[0]
    error = h[1]# yerr argument

    # Taking into account "normed" parameter: 
    if normed == True:
      N = self.G.num_vertices()
      x = [i/N for i in x]
      y = [i/N for i in y]
      error = [i/N for i in error]

    # Taking into account "errobar" parameter and plotting
    if errorbar == True:
      plt.errorbar(x, y, error, fmt="o", color="royalblue", markersize=4) 
    else:
      plt.plot(x, y, "o", color="royalblue", markersize=4)
    plt.show(block=block)
  
  def one_node_cascade(self, fraction_to_fail, initial_node):
    """
    
    one_node_cascade executes failure cascade simulation with the starting failure point equal to the provided initial node (1). 
    Failure cascade algorithm requires going constantly through the network and checking the nodes's statuses (2). 
    The current state of the node is changed to FAILED if the number of node's neighbours that have FAILED statuses exceeds 
    or is equal to fraction_to_fail number (3). Looping over the network finishes when no new FAILED status has been introduced 
    during the iteration (4). The output of the function is the number of nodes with the FAILED status at the end of the simulation (5).

    """

    # Initializing a vector that represents statuses:
    gprop = self.G.new_vertex_property("bool")
    gprop[initial_node] = True #1
    go_on=True

    while go_on == True: #2
      
      go_on=False #4 assume no new FAILED status in the upcoming iteration
      for v in self.G.get_vertices(): #2
        if gprop[v] == 0: # check current node status
          failures = gprop.a[self.G.get_all_neighbors(v)] # extract statuses of all the node's neighbours
          if sum(failures)/len(failures) >= fraction_to_fail:
            gprop[v]=1 #3
            go_on=True # have had new FAILED status, must continue looping

    cascade_size = sum(gprop.a)/len(gprop.a) #5
    
    return (initial_node, cascade_size) 

  def cascade_all_nodes(self, fraction_to_fail = 0.25):
    """

    cascade_all_nodes runs failure cascade simulation (one_node_cascade) for each of the network's nodes to evaluate distribution
    of the final cascade sizes. It returns a dictionary in which each node is assigned a value of the cascade size that it generated.

    """

    nodes_numbers = []
    cascade_sizes =[]

    for v in self.G.get_vertices(): # Take each node
      i, c = self.one_node_cascade(fraction_to_fail, v) # Run for it failure cascade
      nodes_numbers.append(v)#self.G.properties[("v","_graphml_vertex_id")][v]) 
      cascade_sizes.append(c)

    zip_iterator = zip(nodes_numbers, cascade_sizes) # Get pairs of elements.
    dictionary_names_cascade = dict(zip_iterator) # Return dicitionary node_number:cascade_size

    return dictionary_names_cascade 
  
  def plot_cascade(self, dictionary_names_cascade, fraction_to_fail):
    """
    
    plot_cascade generates a histogram for the results of the cascade_all_nodes function.
    It shows the distribution of the failure cascade sizes in the network.
    """

    # Plot settings:
    FONT = 15
    plt.style.use('seaborn-whitegrid')
    plt.figure(figsize=(5,5))
    plt.title("Cascade size histogram C="+ str(fraction_to_fail), fontsize= FONT)
    plt.xlabel("Value", fontsize= FONT)
    plt.ylabel("Fraction of nodes", fontsize= FONT)
               
    # Data transformation for the histogram:
    cascade_sizes = list(dictionary_names_cascade.values())
    unique, counts = np.unique(cascade_sizes, return_counts=True)
    cascade_sizes_counts = dict(zip(unique, counts))
    possible_cascade_sizes, counts = zip(*cascade_sizes_counts.items())
    fractions = [i/sum(counts) for i in counts]

    # Plotting:
    plt.plot(possible_cascade_sizes, fractions,"*", color="royalblue",markersize=4)
    plt.show(block=True)
  
  def robustness_evaluation(self, map_G, step = 1):
    """

    robustness_evaluation performs the robustness measurements according to the provided map_G. 
    Robustness measurements are performed by sorting the nodes according to the map_G values (1).
    Then subsequent fractions of the nodes are taken according to the sorted pattern (2) and removed from the network
    using the filtering option in graph-tool (3). In such a way new subgraphs that contain only not filtered-out (removed) nodes and edges between them
    are generated (4). The largest component sizes of such subnetworks are calculated and returned.

    """

    largest_N = self.G.num_vertices()
    largest_E = self.G.num_edges()
    giant_component_size = []
    vertices_to_remove = map_G.a.argsort()[::-1] # 1
    f_previous = 0
    # settings for a vector that represents whether a node should be taken or not when performing network filtering
    gprop = self.G.new_vertex_property("bool") 
    self.G.vertex_properties["no_removal"] = gprop
    for v in self.G.vertices():
      self.G.properties[("v","no_removal")][v] = True

    for fraction in range(0,100,step):
      f = fraction/100
      new_to_remove = vertices_to_remove[int(f_previous*largest_N):int(f*largest_N)] # (2) adding new nodes to be filtered 

      """ In order to reduce computational costs the filtering statuses are added subsequently. In other words in the first iteration 
      x nodes, equal to f_previous*largest_N, should be filtered (removed), so x nodes have no_removal = False. In new iteration x+y (int(f*largest_N))
      nodes should be added the filtered status. However, already x nodes have no_removal = False, therefore only for node from the range
      int(f_previous*largest_N):int(f*largest_N) must change no_removal = False.
      """

      for node in new_to_remove:
        self.G.properties[("v","no_removal")][node] = False # 3
      
      f_previous = f
      sub = GraphView(self.G, gprop) # 4
      comp, hist = graph_tool.topology.label_components(sub) #5
      giant_component_size.append(max(hist))

    return giant_component_size #5
  
  def plot_robustness(self, metrics_results, step = 1, block = False):
    """

    plot_robustness generates the plots for the data generated by the robustness_evaluation function.

    """

    # Plot settings:
    FONT = 15
    fraction = [i/100 for i in range(0,100,step)]
    plt.figure(figsize = (5,5))
    plt.style.use('seaborn-whitegrid')
    plot_metric_labels = {"Degree": ["-", "#D81B60"] , "Betweenness centrality": ["-", "#1E88E5"],
                          "Closeness centrality" : ["-","#FFC107"], 
                          "Eigenvector centrality": ["-", "#004D40"]} 
    plt.xlabel("Fraction of nodes removed", fontsize = FONT)
    plt.ylabel("Largest component size", fontsize = FONT)
    plt.title("Robustness of the network", fontsize = FONT)

    #Plotting:
    for i in metrics_results:
      data, metric_name = i
      data = [i/max(data) for i in data]
      plt.plot(fraction, data, plot_metric_labels[metric_name][0], label= metric_name, color=plot_metric_labels[metric_name][1], linewidth = 3)

    plt.legend()
    plt.show(block=False)

  def powerlaw(self, cutoff = False):
    """

    powerlaw function adjust the power law distribution according to the Maximum likelihood method for the network's degree sequence.
    The calculations are performed with the usage of poweRlaw library from R package and as the output the value of the adjusted
    alpha parameter is returned. The adjustment is performed for all values from the degree sequence that are larger or equal 
    the cutoff value. If cutoff == False then the cutoff is adjsuted automatically by optimizing the Kolomogrov distance
    between the fitted power law and the data.

    """

    robjects.r('''
        powerlaws <- function(degrees, cutoff = FALSE){
        degrees = as.integer(degrees)
        #print(degrees)

        # Set powerlaw object
        my_powerlaw = displ$new(degrees)

        # Estimate alpha value
        est = estimate_pars(my_powerlaw) 

        # Estimate cutoff value as the one that minimizes Kolomogrov distance between the data and distribution model
        if (cutoff == FALSE){
          est2 = estimate_xmin(my_powerlaw)
          my_powerlaw$setXmin(est2)
          est = estimate_pars(my_powerlaw)
          my_powerlaw$setPars(est$pars)
        }
        else{
          my_powerlaw$setXmin(cutoff)
          est = estimate_pars(my_powerlaw)
          my_powerlaw$setPars(est$pars)
        }

        # Calculate likelihood of the model
        likelihood = dist_ll(my_powerlaw)

        # Calculate percentage of data covered by the powerlaw
        percentage = length(degrees[which(degrees>=my_powerlaw$xmin)])/length(degrees)
        #print(degrees[which(degrees>=my_powerlaw$xmin)])
        
        # Data for plotting the results
        data = plot(my_powerlaw)
        fit = lines(my_powerlaw)
        return(list(data, fit, my_powerlaw$xmin, my_powerlaw$pars, percentage, likelihood, my_powerlaw)) 
        #return(c(my_powerlaw$xmin, my_powerlaw$pars))
        #statistical_test = bootstrap_p(m, no_of_sims = 1000, threads = 2)
        #p_value = statistical_test$p
        
        }''')
    # Make R funtion available for python:
    powerlaw = robjects.globalenv['powerlaws'] 

    # Prepare the degree sequence: 
    degree_map = self.create_degree_distribution_map().fa
    degree_map = degree_map.tolist()

    # Perform calculations:
    power_law_result = powerlaw(degree_map, cutoff) 
    plotting_data = (power_law_result[0][0], power_law_result[0][1], power_law_result[1][0], power_law_result[1][1])
    kmin = power_law_result[2][0]
    alpha = power_law_result[3][0]
    percentage = power_law_result[4][0]
    likelihood = power_law_result[5][0]
    my_powerlaw = power_law_result[6]
    return (kmin, alpha, percentage, likelihood, plotting_data, my_powerlaw)
  
  def bootstrap_powerlaw(self, my_powerlaw, N=100):
    """

    bootstrap_powerlaw calculates the p-value for H0: degree sequence comes from the power law distirbution with parameters: estimated alpha and cutoff; 
    H1: It does not come. The test is performed according to bootstrap_p function from poweRlaw package that simulates N times the data from the distirbution
    and calculates how many times the distance between the theoretical and simulational distributions was larger or equal to the one for the degree sequence.

    """

    robjects.r('''
        assess_p_value <- function(my_powerlaw, N){
        statistical_test = bootstrap_p(my_powerlaw, no_of_sims = N, threads = 2)
        return(statistical_test$p)
        }''')
    p_value = robjects.globalenv['assess_p_value']
    p = p_value(my_powerlaw, N)[0]
    return p
  
  def plot_powerlaw(self, plotting_data, block = False):
    """
    
    plot_powerlaw function visualises the power law fit and the data on the log log scale.

    """
    
    FONT = 15
    # Data preparation:
    datax = plotting_data[0]
    datay = plotting_data[1]
    fitx = plotting_data[2]
    fity = plotting_data[3]
    # Plot settings:
    plt.figure(figsize =(5,5))
    plt.style.use('seaborn-whitegrid')  
    plt.xlabel("log k", fontsize = FONT)
    plt.ylabel("log P(X<k)", fontsize = FONT)
    plt.title("Power law fit", fontsize = FONT)

    # Plotting:
    plt.plot(np.log(datax), np.log(datay), "o", markersize=4, color="#1E88E5")
    plt.plot(np.log(fitx), np.log(fity), linewidth = 3, color = "#FFC107")
    
    plt.show(block = block)
