Module 4: Data Analysis
=======================

Once you have output data from a simulation, we need to be able to *do* something with it. This section will discuss some of the basics of working with simulation data, including some examples from real-world studies. We will discuss various methods of data visualization, as well as the basics of how to apply statistical models to your data. Lastly, we will discuss some best practices for preparing data for publication and sharing so that others can interpret (and reproduce) your simulations as easily as possible.

4.1: Data Visualization
-----------------------

The first step in working with data is to visualize the output so that you can assess system behavior over time (or some other variable of choice). This section will walk through several examples of how to use basic Python scripts to visualize a data set in various ways.

Loading data files
~~~~~~~~~~~~~~~~~~

Numerical data can be loaded from a data file using the ``loadtxt``
function of ``numpy``; i.e., the command is ``np.loadtxt``. You need to
make sure the file is in the same directory as your notebook, or provide
the full path. The filename (or path plus filename) needs to be between
quotes.

Exercise 4.1.#, Loading data and adding a legend
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You are provided with the data files containing the mean montly
temperature of Holland, New York City, and Beijing. The Dutch data is
stored in ``holland_temperature.dat``, and the other filenames are
similar. Plot the temperature for each location against the number of
the month (starting with 1 for January) all in a single graph. Add a
legend by using the function ``plt.legend(['line1','line2'])``, etc.,
but then with more descriptive names. Find out about the ``legend``
command using ``plt.legend?``. Place the legend in an appropriate spot
(the upper left-hand corner may be nice, or let Python figure out the
best place).

.. code:: ipython3

                 
    ! git clone https://github.com/akmadamanchi/ThermoData.git
    
    ### if you get the error "fatal: destination path 'ThermoData' already exists and is not an empty directory."
    ### you can handle this by 1) opening up the menu on the left side of the screen to bring up the table of cotents. 
    ### 2) chose the Files tab in Table of contents.  3) NOTE THIS IS NOT THE File menu at the top of the screen. 
    ### 4) see if there is a folder named ThermoData. 
    ### If there is you can uncomment and run the 'rm -rf ThermoData/' command in the following cell
    

.. code:: ipython3

    #rm -rf ThermoData/ 

.. code:: ipython3

    holland = np.loadtxt('/content/ThermoData/holland_temperature.dat')
    newyork= np.loadtxt('/content/ThermoData/newyork_temperature.dat')
    beijing = np.loadtxt('/content/ThermoData/beijing_temperature.dat')
    plt.plot(np.linspace(1, 12, 12), holland)
    plt.plot(np.linspace(1, 12, 12), newyork)
    plt.plot(np.linspace(1, 12, 12), beijing)
    plt.xlabel('Number of the month')
    plt.ylabel('Mean monthly temperature (Celcius)')
    plt.xticks(np.linspace(1, 12, 12))
    plt.legend(['Holland','New York','Beijing'], loc='best');

Exercise 4.1.#, Subplots and fancy tick markers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Load the average monthly air temperature and seawater temperature for
Holland. Create one plot with two graphs above each other using the
subplot command (use ``plt.subplot?`` to find out how). On the top
graph, plot the air and sea temperature. Label the ticks on the
horizontal axis as ‘jan’, ‘feb’, ‘mar’, etc., rather than 0,1,2,etc. Use
``plt.xticks?`` to find out how. In the bottom graph, plot the
difference between the air and seawater temperature. Add legends, axes
labels, the whole shebang.

Colors
~~~~~~

If you don’t specify a color for a plotting statement, ``matplotlib``
will use its default colors. The first three default colors are special
shades of blue, orange and green. The names of the default colors are a
capital ``C`` followed by the number, starting with number ``0``. For
example

.. code:: ipython3

    plt.plot([0, 1], [0, 1], 'C0')
    plt.plot([0, 1], [1, 2], 'C1')
    plt.plot([0, 1], [2, 3], 'C2')
    plt.legend(['default blue', 'default orange', 'default green']);

There are five different ways to specify your own colors in matplotlib
plotting; you may read about them
`here <http://matplotlib.org/examples/pylab_examples/color_demo.html>`__.
A useful way is to use the html color names. The html codes may be
found, for example, `here <http://en.wikipedia.org/wiki/Web_colors>`__.

.. code:: ipython3

    color1 = 'fuchsia'
    color2 = 'lime'
    color3 = 'DodgerBlue'
    plt.plot([0, 1], [0, 1], color1)
    plt.plot([0, 1], [1, 2], color2)
    plt.plot([0, 1], [2, 3], color3)
    plt.legend([color1, color2, color3]);

The coolest (and nerdiest) way is probably to use the xkcd names, which
need to be prefaced by the ``xkcd:``. The xkcd list of color names is
given by `xkcd <https://xkcd.com/color/rgb/>`__ and includes favorites
such as ‘baby puke green’ and a number of brown colors vary from ``poo``
to ``poop brown`` and ``baby poop brown``. Try it out:

.. code:: ipython3

    plt.plot([1, 2, 3], [4, 5, 2], 'xkcd:baby puke green');
    plt.title('xkcd color baby puke green');

Gallery of graphs
~~~~~~~~~~~~~~~~~

The plotting package ``matplotlib`` allows you to make very fancy
graphs. Check out the matplotlib gallery to get an overview of many of
the options. The following exercises use several of the matplotlib
options.

Exercise 4.1.#, Pie Chart
~~~~~~~~~~~~~~~~~~~~~~~~~

At the 2012 London Olympics, the top ten countries (plus the rest)
receiving gold medals were
``['USA', 'CHN', 'GBR', 'RUS', 'KOR', 'GER', 'FRA', 'ITA', 'HUN', 'AUS', 'OTHER']``.
They received ``[46, 38, 29, 24, 13, 11, 11, 8, 8, 7, 107]`` gold
medals, respectively. Make a pie chart (use ``plt.pie?`` or go to the
pie charts in the matplotlib gallery) of the top 10 gold medal winners
plus the others at the London Olympics. Try some of the keyword
arguments to make the plot look nice. You may want to give the command
``plt.axis('equal')`` to make the scales along the horizontal and
vertical axes equal so that the pie actually looks like a circle rather
than an ellipse. Use the ``colors`` keyword in your pie chart to specify
a sequence of colors. The sequence must be between square brackets, each
color must be between quotes preserving upper and lower cases, and they
must be separated by comma’s like
``['MediumBlue','SpringGreen','BlueViolet']``; the sequence is repeated
if it is not long enough.

Exercise 4.1.#, Fill between
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Load the air and sea temperature, as used in Exercise 4, but this time
make one plot of temperature vs the number of the month and use the
``plt.fill_between`` command to fill the space between the curve and the
horizontal axis. Specify the ``alpha`` keyword, which defines the
transparancy. Some experimentation will give you a good value for alpha
(stay between 0 and 1). Note that you need to specify the color using
the ``color`` keyword argument.

4.2: Statistical Analysis Methods (?)
-------------------------------------

[In Progress: Update Scheduled 06/06/25]

4.3: Model Fitting & Tuning: Examples
-------------------------------------

[In Progress: Update Scheduled 06/09/25]

4.4: Preparing Data for Publication & Sharing
---------------------------------------------

[In Progress: Update Scheduled 06/06/25]