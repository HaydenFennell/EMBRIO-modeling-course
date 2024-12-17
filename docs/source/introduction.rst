Introduction
============

Welcome to the EMBRIO Multiscale Modeling Course! This website will guide you through the fundamentals of creating your own computational biological models, from conceptualizing a biological problem in terms that both you and a computer can understand to performing data fitting and optimization on output from a custom simulation. Computational modeling can be a very complicated subject that often requires knowledge from a variety of different fields. Here, we hope to make the learning process easier by providing a step-by-step introduction to the biological process to help narrow down some of that complexity. No prerequisite background knowledge or experience is required to complete this course; we will cover everything that you need to know to get started. There is much more out there beyond what we provide here, but by the end of this course, you should have the foundational tools that you need to start building models of your specific topics of interest.

Why Learn Modeling?
-------------------

Computational modeling is becoming ubiquitous across virtually all STEM fields. There are many different conceptualizations of modeling across fields, but the core idea is basically always the same. A model is a (simplified) representation of a phenomenon or system that can be used to make predictions (outputs) about the system given certain starting parameters (inputs). Models do not require a computational component, but many of them in modern settings do leverage computational applications. Advancements in computer simulation software and methods over the last several decades has lead to the ability to simulate a broad range of topics and systems that were previously inaccessible. To list a few examples, we can use virtual simulations of biological systems to:

* Understand basic biology of development, homeostasis, and developmental diseases
* Predict effects of chemicals, drugs,… (toxicity) for regulatory applications
* Design new therapies
* Personalize therapies to be more effective for individual patients
* Study systems in which when cells move or tissues reorganize (can be difficult in traditional experimental/imaging contexts)

While we can also use experiments and wet lab work to accomplish many of the same goals, there are a number of benefits to using computational approaches. For example:

* Function as ‘variable magnification’ virtual microscope to allow study of processes at variable spatial and temporal resolution
* Simulations with graded perturbations can infer thresholds for systemic effects when feedback within and between scales can lead to amplification and permanent disruption or compensation and recovery from molecular-scale perturbations
* Infer hard-to-observe parameters and test mechanistic hypotheses
* Reduces need for animal experiments
* Allows testing of situations impractical to test in experiments
* Allows exploration of more combinatorics than experiments 
* Visualize and communicate systems-level outcomes
* A valuable platform for knowledge integration and reuse and to extract more understanding from existing data
* Explore contrafactual conditions not possible in experiment
* Building models codifies and tests biological understanding and reveals critical gaps in understanding and experimental data
* Exploring new conditions and adding replicas is effectively free once you have the models built

If you are not already sold on adding some computational skills to your biological tool belt, we hope to make some of these benefits even more clear over the course of the next several modules. But first it is important that we describe what you can expect to learn.

Learning Objectives
-------------------

Biological modeling is a complicated topic, combining knowledge and expertise from a number of different fields. This course is designed to teach you the basics from each field that you will need to know to get started in building your own biological simulations. By the end of this course, you will be able to:

#. Derive abstract representations of complex biological systems. 
#. Read and interpret basic Python code used in scientific simulations.
#. Create custom Python programs to perform tasks necessary to your system.
#. Navigate the basics of multiscale biological modeling, including (but not limited to):
   
   #. Creating simulations of multicellular phenomenon
   #. Implementing subcellular chemical signalling network models into a multicellular model
   #. Representing the behavior of chemical diffusion and stochastic nutrient delivery in the extracellular matrix.
   
#. Analyze simulation output data to determine fit and validity.
#. Iterate on the model-building process to improve and optimze your models.

How to Use This Site
--------------------

This course will introduce you to the modeling process through 10 Modules of primary content. The modules and a quick summary of their content is provided in the table below.

.. list-table:: 
   :widths: 50 50
   :header-rows: 1
   
   * - Module or Page
     - Description
   * - Module 0: Introduction to Scientific Computing
     - Python basics, loops, data handling
   * - Module 1: Introduction to Biological Modeling
     - The Modeling Workflow, Identifying meaningful hypothesis, model abstraction
   * - Module 2: Building Mathematical Models
     - Modeling with differential equations (ODE/PDE modeling, numerical methods and finite difference solution approaches, common equation types
   * - Module 3: Converting Math To Code
     - Executing equations in code, model specification and parameterization, common simulation methodologies
   * - Module 4: Data Analysis
     - Managing and cleaning data, visualization, model fitting, statistical methods
   * - Module 5: Validity, Reliability, and Optimization
     - Validation vs verification, methods, pairing models with experiment, model optimization
   * - Module 6: Best Practices
     - A collection of best practices, tips, and FAQs about getting started in modeling (suggestions for structuring projects, common pitfalls to avoid, debugging strategies, etc.)
   * - Module 7: Biological background
     - Biological background information needed to complete the worked examples, summary exercises, and projects
   * - Module 8: Worked Examples
     - 
   * - Module 9: Summary Exercises
     - 
   * - Module 10: Projects
     - 

Each module will contain a variety of exercises and coding practice opportunities. Instructions on how to set up a Python coding environment on your computer are included in our manual (see appendix). If you would prefer to use a browser-based implementation of Python, each page is also offered as a fully interactive Jupyter notebook, where you can enter and run your code in cells alongside the instructional material. Each section that offers a Jupyter notebook supplement will have the notebook linked at the top of the section. Please see our appendix on how to use Google Colab for instructions on how to create an editable version of each notebook in your local account.