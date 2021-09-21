import warnings
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
import pandas as pd

class GUI_for_graph_analysis:
  def __init__(self):
    # Initializing the variables and the GUI elements:
    self.G = None
    self.initial_info = widgets.HTML(value = "<b><font color='#555555';font size =5px;font family='Helvetica'>Graphical User Interface for networks analysis</b>")
    self.instruction_header = widgets.HTML(value = "<b><font color='#555555';font size =4px;font family='Helvetica'>Instruction:</b>")
    self.instruction = widgets.HTML(value = "<b><font color='#555555';font size =2.5px;font family='Helvetica'>1. Provide file name with with .graphml extension. <br>2. Hit Prepare the network button. <br>3. Choose the tab of interest. <br>4. Adjust method settings if present.<br> 5. Run the method by hitting the tab's Run button.<br>6. If you want to run new analysis for a new network hit Restart GUI button. </b>")
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
    #self.bootstrap_settings_label = widgets.HTML(value = "<b><font color='black';font size =2px;font family='Helvetica'>Number of simulations for bootstrap </b>") 
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
    #self.cutoff_label = widgets.HTML(value = "<b><font color='black';font size =2px;font family='Helvetica'>Cutoff value</b>")
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
    """self.data_preview_button = widgets.Button(value=False,
                                         description='data_preview_button',
                                         disabled=False,
                                         button_style='', # 'success', 'info', 'warning', 'danger' or ''
                                         tooltip='Description',
                                         icon='check', # (FontAwesome names without the `fa-` prefix)
                                         layout=Layout(width='100%', height='100%'),
                                         style= {'button_color':'#FFD3B4'}
                                         )
    self.data_preview_button.layout.visibility = 'hidden'
    self.data_preview = widgets.Output()"""
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

  def centrality_button_click(self, b):
    """

    Binds the centrality measure choice from the centrality tab with the appropriate map (1) and plot generation (2).

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
          fig, ax = self.G.plot_map_histogram(my_map, self.centrality_choice.value) # 2)
          self.retrieve_data(my_map, "Centrality and clustering")
          my_map = list(my_map.fa)
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

    Defines what to do when the ANND button is clicked. Tutaj dopisacv potem. 

    """
    self.clear()
    if self.error() == True:
      return None
    else:
      corr_value = round(self.G.calculate_assortativity_value(),3)
      corr_meaning = "assortative" if corr_value>0 else "disassortative"
      self.label_corr_value.value = "Degree correlation coefficient equals " + str(corr_value)+". Graph has "+ corr_meaning +' mixing patterns with regards to the degree.' 
      with self.assortativity_out:
        self.assortativity_out.clear_output()
        self.G.plot_ANND(normed = self.ANND_plot_settings_normed.value, errorbar = self.ANND_plot_settings_errorbar.value, block = False)

  def hubs_impact_choice_plot(self, b):
    """
    Binds...
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
            Ns, Es, degrees_set = self.G.hubs_impact_check()
            self.G.plot_hubs_impact1(degrees_set, Es, block = False)

          if self.hubs_impact_choice.value == "Hubs impact 2":
            Ns, Es, degrees_set = self.G.hubs_impact_check()
            self.G.plot_hubs_impact2(degrees_set, Es, Ns, block = False)
  
  def cascade_button_click(self, b):
    self.clear()
    if self.error() == True:
          return None
    else:
      self.button_cascade.style.button_color = '#FFAAA7'
      self.button_cascade.description = "Running simulation..."
      cascade_data = self.G.cascade_all_nodes(fraction_to_fail = self.cascade_fraction_to_fail.value)
      self.retrieve_data(cascade_data, "Cascade")
      with self.cascade_out:
        self.cascade_out.clear_output()
        self.G.plot_cascade(cascade_data, fraction_to_fail = self.cascade_fraction_to_fail.value)
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
          results = self.G.robustness_evaluation(map_G)
          results_to_plot.append([results, metric_name])

      self.retrieve_data(results_to_plot, "Robustness")

      with self.robustness_out:
        self.robustness_out.clear_output()
        self.G.plot_robustness(results_to_plot, block=True)
  
  def powerlaw_button_click(self, b):
    self.clear()
    if self.error() == True:
      return None
    else:
      pvalue = "Not calculated"
      cutoff = self.cutoff.value if self.cutoff_settings.value == False else False
      (kmin, alpha, percentage, likelihood, plotting_data, my_powerlaw) = self.G.powerlaw(cutoff)
      if self.powerlaw_pvalue.value == True:
        # calculate also p-value
        N = self.bootstrap_settings.value
        pvalue = self.G.bootstrap_powerlaw(my_powerlaw, N)
        pvalue = str(round(pvalue, 4))
        self.pvalue_label.layout.visibility = 'visible'
        self.pvalue_value.layout.visibility = 'visible'          

      with self.powerlaw_out:
        self.powerlaw_out.clear_output()
        self.G.plot_powerlaw(plotting_data, block = True)
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
    if self.powerlaw_pvalue.value == True:
      self.bootstrap_settings.layout.visibility = 'visible'
      self.bootstrap_settings_label.layout.visibility = "visible"
    else: 
      self.bootstrap_settings.layout.visibility = 'hidden'
      self.bootstrap_settings_label.layout.visibility = "hidden"
  
  def powerlaw_cutoff(self, b):
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
    display(self.initial_info)
    display(self.instruction_header)
    display(self.instruction)
    preparation = VBox(children = [self.file_name_textbox, self.button_graph_preparation], layout = Layout(width = "100%"))
    display(preparation)
    tabs_preparation = self.tabs
    outs = VBox(children = [self.centrality_out, self.hubs_impact_out, 
                           self.assortativity_out, self.label_corr_value, 
                           self.robustness_out, self.cascade_out, self.powerlaw_out,
                           self.download_button
                         ]) # self.clustering_out
    all = HBox(children = [tabs_preparation, outs])
    #HBox(children=[self.download_button, self.data_preview_button]], layout = Layout(width = "100%")))
    #display(self.file_name_textbox)
    #display(self.button_graph_preparation)
    #display(HBox(children=[self.tabs,
    #                       VBox(children = [
    #                       self.centrality_out, self.hubs_impact_out, 
    #                       self.assortativity_out, self.label_corr_value, 
    #                       self.robustness_out, self.cascade_out]) # self.clustering_out
    #                       ])) # tu na poczatku VBox bylo self.plot_label
    display(all)
    display(self.error_info)
    display(self.restart_button) 
    #display(self.out)

  def bind(self):
    """

    Binds button and other interactivities with the corresponding action functions.

    """

    # Bind prepare graph button with the preparation function:
    self.button_graph_preparation.on_click(self.button_graph_preparation_click)

    # Bind centrality choice with the plot generation and centrality tab
    #self.centrality_out = interactive_output(self.centrality_choice_plot, {'b':self.centrality_choice})
    self.button_centrality.on_click(self.centrality_button_click)
    self.tab_centrality = VBox(children=[self.label_centrality, self.centrality_choice, self.button_centrality])

    # Bind clusterization button with the plot generation and clusterization tab
    #self.button_clustering.on_click(self.clustering_button_click)
    #self.tab_clustering = VBox(children=[self.label_clustering, self.button_clustering])

    # Bind hubs_impact method choice with the plot generation and hubs_impact tab
    #self.hubs_impact_out = interactive_output(self.hubs_impact_choice_plot, {'b':self.hubs_impact_choice})
    self.hubs_impact_button.on_click(self.hubs_impact_choice_plot)
    self.tab_hubs_impact = VBox(children=[self.label_hubs_impact, self.label_hubs_impact_explain, self.hubs_impact_choice, self.hubs_impact_button])

    # Bind assortativity button with the plot generation and assortativity calculation and assortativity tab
    self.button_assortativity.on_click(self.assortativity_button_click)
    self.tab_assortativity = VBox(children=[self.label_ANND_plot, self.label_ANND_plot_settings, 
                                            self.ANND_plot_settings_errorbar, self.ANND_plot_settings_normed, self.button_assortativity
                                            ])

    self.button_robustness.on_click(self.robustness_button_click)
    self.robustness = VBox(children=[self.label_robustness_info, self.label_robustness_settings, self.robustness_degree, self.robustness_betweenness, 
                                     self.robustness_closeness,
                                     self.robustness_eigenvector, self.button_robustness])

    self.button_cascade.on_click(self.cascade_button_click)
    self.tab_cascade = VBox(children=[self.cascade_info, HBox(children = [self.cascade_fraction_to_fail_label, self.cascade_fraction_to_fail]), 
                                      self.button_cascade])
    
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

    #self.data_preview_button.on_click(self.data_preview_button_click)
    self.restart_button.on_click(self.gui_restart)

  def gui_restart(self,b):
    self.G = None
    self.file_name_textbox.value = "Provide file name here"
    self.button_graph_preparation.description = "Prepare the graph"
    self.button_graph_preparation.style.button_color = "#FFAAA7"
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
    # zmienic dac do funkcji on click
    if self.G == None or self.file_name_textbox.value == "No file name provided. Provide file name here." or self.file_name_textbox.value == "":
      if return_message==True:
        self.error_info.value = "<b><font color='#FFAAA7';font size =3px;font family='Helvetica'>Cannot use the method. Provide file name and prepare the network first.</b>"
      return True
  
  def clear(self):
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
      self.dataframe = pd.DataFrame({"NodeIndex":nodes, "Cascade size": values})
      self.download_button.layout.visibility = 'visible'
      self.download_button.contents = lambda: self.dataframe.to_string()
  
  #def data_preview_button_click(self, b):
   # with self.data_preview:
    #  print(self.dataframe.head(n=5))
  

