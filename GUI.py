import warnings
warnings.filterwarnings("ignore")
from graph_tool.all import *
import graph_tool.all as gt
import ipywidgets as widgets
import matplotlib.pyplot as plt
import numpy as np
from IPython.display import clear_output
import seaborn as sns
from ipywidgets import *
import rpy2.robjects.packages as rpackages
from rpy2.robjects.packages import importr
import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.vectors import StrVector
import pandas as pd
from IPython.utils import io

# Libraries for Download button
import base64
import hashlib
from typing import Callable

import ipywidgets
from IPython.display import HTML, display

# Installing R packages
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

# Defining additional ipywidget that will perform data download after button hitting - DownloadButton
class DownloadButton(ipywidgets.Button):
    """Download button with dynamic content

    The content is generated using a callback when the button is clicked. It is defined as an extension of "button" class in ipywidgets (source: https://stackoverflow.com/questions/61708701/how-to-download-a-file-using-ipywidget-button). 
    """

    def __init__(self, filename: str, contents: Callable[[], str], **kwargs):
        super(DownloadButton, self).__init__(**kwargs)
        self.filename = filename
        self.contents = contents
        self.on_click(self.__on_click)

    def __on_click(self, b):
        contents: bytes = self.contents().encode('utf-8')
        b64 = base64.b64encode(contents)
        payload = b64.decode()
        digest = hashlib.md5(contents).hexdigest()  # bypass browser cache
        id = f'dl_{digest}'

        display(HTML(f"""
<html>
<body>
<a id="{id}" download="{self.filename}" href="data:text/csv;base64,{payload}" download>
</a>

<script>
(function download() {{
document.getElementById('{id}').click();
}})()
</script>

</body>
</html>
"""))

# Graphical User Interface: 
class GUI_for_network_analysis:
  def __init__(self):
    # Initializing the variables and the GUI elements:
    self.G = None
    self.initial_info = widgets.HTML(value = "<b><font color='#555555';font size =5px;font family='Helvetica'>Graphical User Interface for networks analysis</b>")
    self.instruction_header = widgets.HTML(value = "<b><font color='#555555';font size =4px;font family='Helvetica'>Instruction:</b>")
    self.instruction = widgets.HTML(value = "<b><font color='#555555';font size =2.5px;font family='Helvetica'>1. Provide file name with with .graphml extension. <br>2. Hit Prepare the network button (Parallel links, nodes not from the largest component will be removed. Network is set as undirected) . <br>3. Choose the tab of interest. <br>4. Adjust method settings if present.<br> 5. Run the method by hitting the tab's Run button.<br>6. If you want to run new analysis for a new network hit Restart GUI button. </b>")
    self.file_name_textbox = widgets.Text(value='Provide file name here',
                                          placeholder='Type something',
                                          description='Network:',
                                          disabled=False,
                                          align_items='center',
                                          layout=Layout(width='40%')#, height='10px')
                                          )
    self.button_graph_preparation = widgets.Button(value=False,
                                                   description='Prepare the network',
                                                   disabled=False,
                                                   button_style='', # 'success', 'info', 'warning', 'danger' or ''
                                                   tooltip='Description',
                                                   icon='check', # (FontAwesome names without the `fa-` prefix)
                                                   layout=Layout(width='40%', height='20%'),
                                                   style= {'button_color':'#FFAAA7'}
                                                   )
    self.links_nodes_number_info = widgets.Label(value="")
      
    self.label_centrality = widgets.HTML(value = "<b><font color='black';font size =2px;font family='Helvetica'>Histograms of centrality measures</b>") 
    self.centrality_choice = widgets.Dropdown(
        options=['Choose from the list','Degree', 'Betweenness centrality', 'Closeness centrality',
                 'Eigenvector centrality', "Clustering coefficient"],
        description='Measure: ',
        disabled=False,
        layout=Layout(width='90%')
        )
    self.button_centrality = widgets.Button(value=False,
                                            description='Run',
                                            disabled=False,
                                            button_style='', # 'success', 'info', 'warning', 'danger' or ''
                                            tooltip='Description',
                                            icon='check', # (FontAwesome names without the `fa-` prefix)
                                            layout=Layout(width='90%', height='20%'),
                                            style= {'button_color':'#98DDCA'}
                                            )
    self.centrality_out = widgets.Output()
    self.info_mini = widgets.HTML(value = "<b><font color='black';font size =2px;font family='Helvetica'>Minimum: </b>")
    self.info_mini_value = widgets.Label(value = "")
    self.info_maxi = widgets.HTML(value = "<b><font color='black';font size =2px;font family='Helvetica'>Maximum: </b>")
    self.info_maxi_value = widgets.Label(value = "")
    self.info_avg = widgets.HTML(value = "<b><font color='black';font size =2px;font family='Helvetica'>Average: </b>")
    self.info_avg_value = widgets.Label(value = "")
    self.info_std = widgets.HTML(value="<b><font color='black';font size =2px;font family='Helvetica'>Standard deviation: </b>")
    self.info_std_value = widgets.Label(value = "")
    
    self.button_assortativity = widgets.Button(value=False,
                                                   description='Run',
                                                   disabled=False,
                                                   button_style='', # 'success', 'info', 'warning', 'danger' or ''
                                                   tooltip='Description',
                                                   icon='check', # (FontAwesome names without the `fa-` prefix)
                                                   layout=Layout(width='90%', height='20%'),
                                                   style= {'button_color':'#98DDCA'}
                                                   ) #można zrobić pogrubione (działa) , "font_weight":"bold" dodać do stylu

    self.label_corr_value = widgets.Label(value = "")  #było " "
    self.label_ANND_plot = widgets.HTML(value = "<b><font color='black';font size =2px;font family='Helvetica'>Assortativity examination: Average Nearest Neighbour Degree (ANND) plot and degree correlation coefficient</b>") 
    self.label_ANND_plot_settings = widgets.Label(value = "ANND plot settings:")
    self.ANND_plot_settings_normed = widgets.Checkbox(value=False,
                                                      description='Normed ANND',
                                                      disabled=False,
                                                      indent=False)
    self.ANND_plot_settings_errorbar = widgets.Checkbox(value=False,
                                                      description='Errorbars',
                                                      disabled=False,
                                                      indent=False)
    self.assortativity_out = widgets.Output()

    self.hubs_impact_choice = widgets.Dropdown(
        options=['Choose from the list','Hubs impact 1', 'Hubs impact 2'],
        description='Measure: ',
        disabled=False,
        layout=Layout(width='90%')
        )
    self.hubs_impact_button = widgets.Button(value=False,
                                                   description='Run',
                                                   disabled=False,
                                                   button_style='', # 'success', 'info', 'warning', 'danger' or ''
                                                   tooltip='Description',
                                                   icon='check', # (FontAwesome names without the `fa-` prefix)
                                                   layout=Layout(width='90%', height='20%'),
                                                   style= {'button_color':'#98DDCA'}
                                                   )
    self.label_hubs_impact = widgets.HTML(value = "<b><font color='black';font size =2px;font family='Helvetica'>Examination of hubs impact</b>")
    self.label_hubs_impact_explain = widgets.Label(value = "Hubs impact examination consists of creating subnetworks.. i tutaj walnąć ten ładny matematyczny zapis z mgr")
    self.hubs_impact_out = widgets.Output()

    self.button_robustness = widgets.Button(value=False,
                                                   description='Run',
                                                   disabled=False,
                                                   button_style='', # 'success', 'info', 'warning', 'danger' or ''
                                                   tooltip='Description',
                                                   icon='check', # (FontAwesome names without the `fa-` prefix)
                                                   layout=Layout(width='90%', height='20%'),
                                                   style= {'button_color':'#98DDCA'}
                                                   )
    self.robustness_degree = widgets.Checkbox(value=True,
                                              description='Degree',
                                              disabled=False,
                                              indent=False)
    self.robustness_betweenness = widgets.Checkbox(value=False,
                                              description='Betweennness centrality',
                                              disabled=False,
                                              indent=False)
    self.robustness_closeness = widgets.Checkbox(value=False,
                                              description='Closeness centrality',
                                              disabled=False,
                                              indent=False)
    self.robustness_eigenvector = widgets.Checkbox(value=False,
                                              description='Eigenvector centrality',
                                              disabled=False,
                                              indent=False)
    self.label_robustness_info = widgets.HTML(value = "<b><font color='black';font size =2px;font family='Helvetica'>Examination of the network robustness</b>") 
    self.label_robustness_settings = widgets.Label(value = "Choose metrics for the network robustness examination:") 
    self.robustness_out = widgets.Output()
    
    self.cascade_info = widgets.HTML(value = "<b><font color='black';font size =2px;font family='Helvetica'>Simulation of failure cascade</b>") 
    self.button_cascade = widgets.Button(value=False,
                                         description='Run',
                                         disabled=False,
                                         button_style='', # 'success', 'info', 'warning', 'danger' or ''
                                         tooltip='Description',
                                         icon='check', # (FontAwesome names without the `fa-` prefix)
                                         layout=Layout(width='90%', height='20%'),
                                         style= {'button_color':'#98DDCA'}
                                         )
    self.cascade_fraction_to_fail = widgets.FloatSlider(value=0.25, min=0, max=1, step=0.05, 
                                                      description='',
                                                      disabled=False,
                                                      continuous_update=False,
                                                      orientation='horizontal',
                                                      readout=True,
                                                      readout_format='.2f')
    self.cascade_fraction_to_fail_label = widgets.Label(value = "Failure fraction")
    self.cascade_out = widgets.Output()

    self.label_powerlaw = widgets.HTML(value = "<b><font color='black';font size =2px;font family='Helvetica'>Fitting power law to the degree sequence using Maximum Likelihood estimator</b>")
    self.powerlaw_settings = widgets.HTML(value = "Settings:")
    self.powerlaw_pvalue = widgets.Checkbox(value=False,
                                      description='Calculate p-value',
                                      disabled=False,
                                      indent=False)
    self.bootstrap_settings_label = widgets.Label(value = "Number of simulations for bootstrap")
    self.bootstrap_settings = widgets.IntSlider(value=100, min=10, max=1000, step=100, 
                                                      description='',
                                                      disabled=False,
                                                      continuous_update=False,
                                                      orientation='horizontal',
                                                      readout=True
                                                      )
    self.bootstrap_settings.layout.visibility = 'hidden'
    self.bootstrap_settings_label.layout.visibility = 'hidden'
    self.cutoff_settings = widgets.Checkbox(value=True,
                                      description='Cutoff value according to Kolomogrov distance',
                                      disabled=False,
                                      indent=False)
    self.cutoff_label = widgets.Label(value = "Cutoff value") 
    self.cutoff_label.layout.visibility = 'hidden'
    self.cutoff = widgets.IntSlider(value = 1, min=1, max=100, step=1, 
                                                      description='',
                                                      disabled=False,
                                                      continuous_update=False,
                                                      orientation='horizontal',
                                                      readout=True
                                                      )
    self.cutoff.layout.visibility = 'hidden'
    self.pvalue_label = widgets.HTML(value = "<b><font color='black';font size =2px;font family='Helvetica'>P-value:</b>")
    self.pvalue_value = widgets.Label(value="")
    self.pvalue_label.layout.visibility = 'hidden'
    self.pvalue_value.layout.visibility = 'hidden'
    self.powerlaw_button = widgets.Button(value=False,
                                         description='Run',
                                         disabled=False,
                                         button_style='', # 'success', 'info', 'warning', 'danger' or ''
                                         tooltip='Description',
                                         icon='check', # (FontAwesome names without the `fa-` prefix)
                                         layout=Layout(width='90%', height='20%'),
                                         style= {'button_color':'#98DDCA'}
                                         )
    self.powerlaw_out = widgets.Output()
    
    self.restart_button = widgets.Button(value=False,
                                         description='Restart GUI',
                                         disabled=False,
                                         button_style='', # 'success', 'info', 'warning', 'danger' or ''
                                         tooltip='Description',
                                         icon='check', # (FontAwesome names without the `fa-` prefix)
                                         layout=Layout(width='40%', height='20%'),
                                         style= {'button_color':'#FFD3B4'}
                                         )

    self.error_info = widgets.HTML(value = " ")
    self.plot_label = widgets.HTML(value = "Plot and info")
    self.download_button = DownloadButton(filename='data.csv', contents=lambda: f'', description='Download data')
    self.download_button.layout.visibility = 'hidden'
    self.download_button.layout.width = '90%'
    self.download_button.style.button_color = '#D5ECC2'
    self.dataframe = None 
    

  def button_graph_preparation_click(self, button):
    """

    Defines what to do when the graph preparation button is clicked.

    """
    self.clear()
    
    # Error handling:
    if self.file_name_textbox.value == "" or self.file_name_textbox.value == 'Provide file name here':
      self.file_name_textbox.value = "No file name provided. Provide file name here."
      return None
    if ".graphml" not in self.file_name_textbox.value:
      self.file_name_textbox.value = "Incorrect file name. File must have .graphml extension."
      return None

    self.button_graph_preparation.description = "Preparing..."
    self.error_info.value = " "

    # Graph upload from the file:
    self.G = My_Network(self.file_name_textbox.value)

    # Graph preparation - removal of the parallel edges, non-connected components etc.:
    self.G.prepare_the_network()

    self.button_graph_preparation.description = "Graph is ready! Now choose the tool below."
    self.button_graph_preparation.style.button_color = '#D5ECC2'
    self.links_nodes_number_info.value = "Number of nodes: "+str(self.G.G.num_vertices())+", Number of links: " + str(self.G.G.num_edges())

  def centrality_button_click(self, b):
    """

    Binds the centrality measure button from the centrality tab with the appropriate map (1) and plot generation (2) and statistics calculations (3).

    """
    self.clear()
    with self.centrality_out:
      if self.centrality_choice.value == "Choose from the list":
        pass
      else:
        # 1):
        if self.error() == True:
          return None
        else:
          centrality_choices_functions = {'Degree':self.G.create_degree_distribution_map, 
                                          'Betweenness centrality':self.G.create_betweenness_distribution_map,
                                          'Closeness centrality': self.G.create_closeness_distribution_map, 
                                          'Eigenvector centrality':self.G.create_eigenvector_distribution_map,
                                          "Clustering coefficient": self.G.create_clustering_map}
          my_map = centrality_choices_functions[self.centrality_choice.value]() 
          fig, ax = self.G.plot_map_histogram(my_map, self.centrality_choice.value) # 2
          self.retrieve_data(my_map, "Centrality and clustering")
          my_map = list(my_map.fa)
          # 3:
          self.info_mini_value.value = str(min(my_map))
          self.info_maxi_value.value = str(max(my_map))
          self.info_avg_value.value = str(round(np.mean(my_map),4))
          self.info_std_value.value = str(round(np.std(my_map),4))
          self.info_mini = widgets.HTML(value = "<b><font color='black';font size =2px;font family='Helvetica'>Minimum: </b>")
          self.info_maxi = widgets.HTML(value = "<b><font color='black';font size =2px;font family='Helvetica'>Maximum: </b>")
          self.info_avg = widgets.HTML(value = "<b><font color='black';font size =2px;font family='Helvetica'>Average: </b>")
          self.info_std = widgets.HTML(value="<b><font color='black';font size =2px;font family='Helvetica'>Standard deviation: </b>")
          display(VBox(children = [
              HBox(children= [self.info_mini, self.info_mini_value]),
              HBox(children= [self.info_maxi, self.info_maxi_value]),
              HBox(children= [self.info_avg, self.info_avg_value]),
              HBox(children= [self.info_std, self.info_std_value]) 
          ]))

  
  def assortativity_button_click(self, b):
    """

    Binds the assortativity button with the ANND plot generation (1) and degree correlation calculations (2).

    """
    self.clear()
    if self.error() == True:
      return None
    else:
      corr_value = round(self.G.calculate_assortativity_value(),3)
      corr_meaning = "assortative" if corr_value>0 else "disassortative"
      self.label_corr_value.value = "Degree correlation coefficient equals " + str(corr_value)+". Graph has "+ corr_meaning +' mixing patterns with regards to the degree.' # 2
      with self.assortativity_out:
        self.assortativity_out.clear_output()
        self.G.plot_ANND(normed = self.ANND_plot_settings_normed.value, errorbar = self.ANND_plot_settings_errorbar.value, block = False) # 1

  def hubs_impact_choice_plot(self, b):
    """
    
    Binds the hubs impact button with the hubs impact plot generation. Data is firstly calculated by calling hubs_impact check function (1) and then plotted (2).
    
    """
    self.clear()
    with self.hubs_impact_out:
      if self.hubs_impact_choice.value == "Choose from the list":
        pass
      else:
        if self.error() == True:
          return None
        else:
          if self.hubs_impact_choice.value == "Hubs impact 1":
            Ns, Es, degrees_set = self.G.hubs_impact_check() # 1
            self.G.plot_hubs_impact1(degrees_set, Es, block = False) # 2

          if self.hubs_impact_choice.value == "Hubs impact 2":
            Ns, Es, degrees_set = self.G.hubs_impact_check() # 1
            self.G.plot_hubs_impact2(degrees_set, Es, Ns, block = False) # 2
  
  def cascade_button_click(self, b):
    """
    
    Binds the cascade button with fialure cascade simulation performance (1), plotting (2) and the statistics calculations (3). 
    
    """
    self.clear()
    if self.error() == True:
          return None
    else:
      # Button settings: 
      self.button_cascade.style.button_color = '#FFAAA7'
      self.button_cascade.description = "Running simulation..."
      
      # Data generation: 
      cascade_data = self.G.cascade_all_nodes(fraction_to_fail = self.cascade_fraction_to_fail.value) # 1
      self.retrieve_data(cascade_data, "Cascade")
      with self.cascade_out:
        self.cascade_out.clear_output()
        self.G.plot_cascade(cascade_data, fraction_to_fail = self.cascade_fraction_to_fail.value) # 2
        # (3):
        self.info_mini_value.value = str(min(cascade_data.values()))
        self.info_maxi_value.value = str(max(cascade_data.values()))
        self.info_avg_value.value = str(round(np.mean(list(cascade_data.values())),4))
        self.info_std_value.value = str(round(np.std(list(cascade_data.values())),4))
        self.info_mini = widgets.HTML(value = "<b><font color='black';font size =2px;font family='Helvetica'>Minimum: </b>")
        self.info_maxi = widgets.HTML(value = "<b><font color='black';font size =2px;font family='Helvetica'>Maximum: </b>")
        self.info_avg = widgets.HTML(value = "<b><font color='black';font size =2px;font family='Helvetica'>Average: </b>")
        self.info_std = widgets.HTML(value="<b><font color='black';font size =2px;font family='Helvetica'>Standard deviation: </b>")
        display(VBox(children = [
            HBox(children= [self.info_mini, self.info_mini_value]),
            HBox(children= [self.info_maxi, self.info_maxi_value]),
            HBox(children= [self.info_avg, self.info_avg_value]),
            HBox(children= [self.info_std, self.info_std_value]) 
        ]))
        
      self.button_cascade.description = "Run failure cascade simulation"
      self.button_cascade.style.button_color = '#98DDCA'
  
  def robustness_button_click(self, b):
    """
    
    Binds robustness button with the reboustness examination. In the call the data is generated (1) and then plotted (2).
    
    """
    self.clear()
    if self.error() == True:
          return None
    else:
      metrics_to_run = {self.robustness_degree:[self.G.create_degree_distribution_map, "Degree"], 
                        self.robustness_betweenness:[self.G.create_betweenness_distribution_map, "Betweenness centrality"] , 
                        self.robustness_closeness:[self.G.create_closeness_distribution_map, 'Closeness centrality'],
                        self.robustness_eigenvector:[self.G.create_eigenvector_distribution_map,'Eigenvector centrality']}

      results_to_plot = []
      for metric in metrics_to_run.keys():
        if metric.value == True:
          [function, metric_name] = metrics_to_run[metric]
          map_G = function()
          results = self.G.robustness_evaluation(map_G) # 1
          results_to_plot.append([results, metric_name])

      self.retrieve_data(results_to_plot, "Robustness")

      with self.robustness_out:
        self.robustness_out.clear_output()
        self.G.plot_robustness(results_to_plot, block=True) # 2
  
  def powerlaw_button_click(self, b):
    """
    
    Binds the powerlaw button with the power law adjustment to the degree sequence. Parameters are calculated (1), the fit is plotted (2) and the statistics are calculated (3).
    """
    self.clear()
    if self.error() == True:
      return None
    else:
      pvalue = "Not calculated"
      cutoff = self.cutoff.value if self.cutoff_settings.value == False else False
      (kmin, alpha, percentage, likelihood, plotting_data, my_powerlaw) = self.G.powerlaw(cutoff) # 1
      if self.powerlaw_pvalue.value == True:
        # calculate also p-value
        N = self.bootstrap_settings.value
        pvalue = self.G.bootstrap_powerlaw(my_powerlaw, N)
        pvalue = str(round(pvalue, 4))
        self.pvalue_label.layout.visibility = 'visible'
        self.pvalue_value.layout.visibility = 'visible'          

      with self.powerlaw_out:
        self.powerlaw_out.clear_output()
        self.G.plot_powerlaw(plotting_data, block = True) # 2
        # 3:
        self.info_mini.value = "<b><font color='black';font size =2px;font family='Helvetica'>Cutoff: </b>"
        self.info_mini_value.value = str(kmin)
        self.info_maxi.value = "<b><font color='black';font size =2px;font family='Helvetica'>Power law parameter alpha: </b>"
        self.info_maxi_value.value = str(round(alpha,4))
        if alpha>3 or alpha<2:
          self.info_maxi_value.value+= ", anomalous regime, should be 2<alpha<3"
        self.info_avg.value = "<b><font color='black';font size =2px;font family='Helvetica'>Percentage of data covered: </b>"
        self.info_avg_value.value = str(round(percentage*100,4))
        self.info_std.value = "<b><font color='black';font size =2px;font family='Helvetica'>Likelihood: </b>"
        self.info_std_value.value = str(round(likelihood,4))
        self.pvalue_value.value = pvalue
        display(VBox(children = [
            HBox(children= [self.info_mini, self.info_mini_value]),
            HBox(children= [self.info_maxi, self.info_maxi_value]),
            HBox(children= [self.info_std, self.info_std_value]),
            HBox(children= [self.info_avg, self.info_avg_value]),
            HBox(children= [self.pvalue_label, self.pvalue_value])
        ]))
  
  def powerlaw_pvalue_true(self, b):
    """
    
    Function for handling the powerlaw settings. It makes visible the bootstrap settings if the pvalue is to be assesed (pvalue checkbox is True).
    
    """
    if self.powerlaw_pvalue.value == True:
      self.bootstrap_settings.layout.visibility = 'visible'
      self.bootstrap_settings_label.layout.visibility = "visible"
    else: 
      self.bootstrap_settings.layout.visibility = 'hidden'
      self.bootstrap_settings_label.layout.visibility = "hidden"
  
  def powerlaw_cutoff(self, b):
    """
    
    Function for handling the powerlaw settings. It makes visible the cutoff choice bar if the default option for cutoff adjustment using the Kolomogrov distance is not chosen.
    
    """
    if self.cutoff_settings.value == False:
      self.cutoff_label.layout.visibility = "visible"
      self.cutoff.layout.visibility = 'visible'
      if self.error(return_message = False) == True:
        return None
      else:
        degree_values = self.G.create_degree_distribution_map().fa
        self.cutoff.min = min(degree_values)
        self.cutoff.max = max(degree_values)
        self.cutoff.value = self.cutoff.min
    else:
      self.cutoff_label.layout.visibility = "hidden"
      self.cutoff.layout.visibility = 'hidden'

  def display(self):
    """
    
    Displays all the elements of the GUI in the appropriate order to form the interface.
    
    """
    display(self.initial_info)
    display(self.instruction_header)
    display(self.instruction)
    preparation = VBox(children = [self.file_name_textbox, self.button_graph_preparation, self.links_nodes_number_info], layout = Layout(width = "100%"))
    display(preparation)
    tabs_preparation = self.tabs
    outs = VBox(children = [self.centrality_out, self.hubs_impact_out, 
                           self.assortativity_out, self.label_corr_value, 
                           self.robustness_out, self.cascade_out, self.powerlaw_out,
                           self.download_button
                         ]) # self.clustering_out
    all = HBox(children = [tabs_preparation, outs])
    display(all)
    display(self.error_info)
    display(self.restart_button) 

  def bind(self):
    """

    Binds button and other interactivities with the corresponding action functions.

    """

    # Bind prepare graph button with the preparation function:
    self.button_graph_preparation.on_click(self.button_graph_preparation_click)

    # Bind centrality choice button with the centrality examination and centrality tab
    self.button_centrality.on_click(self.centrality_button_click)
    self.tab_centrality = VBox(children=[self.label_centrality, self.centrality_choice, self.button_centrality])

    # Bind hubs_impact button with the plot generation and hubs_impact tab
    self.hubs_impact_button.on_click(self.hubs_impact_choice_plot)
    self.tab_hubs_impact = VBox(children=[self.label_hubs_impact, self.label_hubs_impact_explain, self.hubs_impact_choice, self.hubs_impact_button])

    # Bind assortativity button with the assortativity examination and assortativity tab
    self.button_assortativity.on_click(self.assortativity_button_click)
    self.tab_assortativity = VBox(children=[self.label_ANND_plot, self.label_ANND_plot_settings, 
                                            self.ANND_plot_settings_errorbar, self.ANND_plot_settings_normed, self.button_assortativity
                                            ])
    
    # Bind robustness button with the robustness examination and robustness tab
    self.button_robustness.on_click(self.robustness_button_click)
    self.robustness = VBox(children=[self.label_robustness_info, self.label_robustness_settings, self.robustness_degree, self.robustness_betweenness, 
                                     self.robustness_closeness,
                                     self.robustness_eigenvector, self.button_robustness])
    
    # Bind cascade button with the failure cascade examination and cascade tab
    self.button_cascade.on_click(self.cascade_button_click)
    self.tab_cascade = VBox(children=[self.cascade_info, HBox(children = [self.cascade_fraction_to_fail_label, self.cascade_fraction_to_fail]), 
                                      self.button_cascade])
    
    # Bind powerlaw button with the powerlaw examination, bind powerlaw settings with the corresponding actions, add all to the powerlaw tab
    self.powerlaw_button.on_click(self.powerlaw_button_click)
    self.powerlaw_bootstrap = interactive_output(self.powerlaw_pvalue_true, {'b':self.powerlaw_pvalue})
    self.powerlaw_cutoff = interactive_output(self.powerlaw_cutoff, {'b':self.cutoff_settings})
    self.tab_powerlaw = VBox(children = [self.label_powerlaw, self.powerlaw_settings, self.powerlaw_pvalue,
                                         self.powerlaw_bootstrap,
                                         self.bootstrap_settings_label, self.bootstrap_settings, 
                                         self.powerlaw_cutoff, self.cutoff_settings, self.cutoff_label, 
                                         self.cutoff,
                                         self.powerlaw_button])

    # Joining tabs in the GUI
    self.tabs = widgets.Accordion(children = [self.tab_centrality, 
                                             self.tab_assortativity, 
                                              self.tab_hubs_impact, self.robustness, self.tab_cascade, self.tab_powerlaw], 
                                  layout=Layout(width='40%', min_width = "300px",
                                                ), selected_index = None) #self.tab_clustering bylo kiedys, 
                                                #layout in_height='500px',max_height='500px',  display='flex'align_items='stretch'
  
    
    # Additional tabs' settings
    self.tabs.set_title(0, '> Centrality and clusterization ')
    #self.tabs.set_title(1, "> Clusterization")
    self.tabs.set_title(1, '> Assortativity')
    self.tabs.set_title(2, '> Hubs impact')
    self.tabs.set_title(3, '> Robustenss')
    self.tabs.set_title(4, '> Failure cascade')
    self.tabs.set_title(5, '> Power law fitting')

    # Bind restart button with the restart function
    self.restart_button.on_click(self.gui_restart)

  def gui_restart(self,b):
    """ 
    
    Sets everything to the initial settings by cleaning the output widgets, fixing colors, bringing original texts to the labels and buttons
    
    """
    self.G = None
    self.file_name_textbox.value = "Provide file name here"
    self.button_graph_preparation.description = "Prepare the graph"
    self.button_graph_preparation.style.button_color = "#FFAAA7"
    self.links_nodes_number_info.value = ""
    self.centrality_choice.value = "Choose from the list"
    self.centrality_out.clear_output()
    #self.clustering_out.clear_output()
    self.hubs_impact_choice.value = "Choose from the list"
    self.hubs_impact_out.clear_output()
    self.label_corr_value.value = ""
    self.ANND_plot_settings_normed.value = False
    self.ANND_plot_settings_errorbar.value = False
    self.assortativity_out.clear_output() 
    self.cascade_fraction_to_fail.value = 0.25
    self.cascade_out.clear_output()
    self.robustness_degree.value = False
    self.robustness_betweenness.value = False
    self.robustness_closeness.value = False
    self.robustness_eigenvector.value = False
    self.robustness_out.clear_output()
    self.powerlaw_pvalue.value = False
    self.cutoff_settings.value = True
    #self.data_preview.clear_output()
    #self.data_preview_button.layout.visibility = 'hidden'
    self.download_button.layout.visibility = 'hidden'
  
  def error(self, return_message = True):
    """
    
    Used for error handling - checks if the file is provided in the appropriate format. This functions is called always before running any of the methods in the GUI.
    
    """
    if self.G == None or self.file_name_textbox.value == "No file name provided. Provide file name here." or self.file_name_textbox.value == "":
      if return_message==True:
        self.error_info.value = "<b><font color='#FFAAA7';font size =3px;font family='Helvetica'>Cannot use the method. Provide file name and prepare the network first.</b>"
      return True
  
  def clear(self):
    """
    
    Clears the outputs. Used to make previous plots and statistics disappear from the GUI when the new method is called. This functions is called always before running any of the methods in the GUI.
    
    """
    self.centrality_out.clear_output() 
    self.hubs_impact_out.clear_output()
    self.assortativity_out.clear_output()
    self.robustness_out.clear_output()
    #self.clustering_out.clear_output()
    self.cascade_out.clear_output()
    self.powerlaw_out.clear_output()
    self.label_corr_value.value = ""
    #self.data_preview.clear_output()
    #self.data_preview_button.layout.visibility = 'hidden'
    self.download_button.layout.visibility = 'hidden'
  
  def retrieve_data(self, data, method):
    """
    
    Used to gather the data from the method functions so that it is downloadable. Called in 3 cases - when the robustness, cascade or Centrality and clustering methods are chosen.
    """
    if method == "Centrality and clustering":
      my_map = data
      my_map_values = my_map.a[self.G.G.get_vertices()]
      nodes = self.G.G.get_vertices()
      self.dataframe = pd.DataFrame({"NodeIndex":nodes, "MeasureValue": my_map_values})
      #self.data_preview_button.layout.visibility = 'visible'
      self.download_button.layout.visibility = 'visible'
      self.download_button.contents = lambda: self.dataframe.to_string()

    if method == "Robustness":
      results_to_plot = data
      dataframe = {}
      for row in results_to_plot:
        dataframe[row[1]] = row[0]
      self.dataframe = pd.DataFrame(dataframe)
      self.download_button.layout.visibility = 'visible'
      self.download_button.contents = lambda: self.dataframe.to_string()

    if method == "Cascade":
      nodes = data.keys()
      values = data.values()
      fractions = [i/100 for i in range(0,100)]
      self.dataframe = pd.DataFrame({"Fraction": fractions, "NodeIndex":nodes, "Cascade size": values})
      self.download_button.layout.visibility = 'visible'
      self.download_button.contents = lambda: self.dataframe.to_string()
  

  

