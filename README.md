Plug-in tool within ABAQUS CAE to perform energy and stength calculations along any given desired tensile fracture plane

There are three main files within the plugin:

1. EDC2.py -
This is the code responsible for performing all calculations to determine the energy and strength of desired tensile fracture planes.
	- Within this file there are multiple functions that are used:
		- axisplanedel: deletes all previous graphics produced by FracTool simulations
		- plycomp: if a hexahedral mesh is chosen this finds the cut plane area in each element split by the number of plies
		- plycomptet: if a tetrahedral mesh is chosen this finds the cut plane area in each element split by the number of plies
		- plycomphextet: if a mixed hexahedral and tetrahedral mesh is chosen this finds the cut plane area in each element split by the number of plies
		- plycomphexel: for a mixed mesh this find the cut plane area in each element for hexahedral elements
		- plycomptetel: for a mixed mesh this find the cut plane area in each element for tetrahedral elements 
		- line: sorts points for the tetrahedral mesh
		- poly4: calculates area of quadrilateral with 3D nodes
		- poly5: calculates area of polygons with more than 4 3D nodes
		- fibremat: calculates failure energy and force for given plane
		- plot: plots graphs and creates Results.csv file with results of running the plug-in

2. plugin3DB.py - 
This code is responsible for creating the user interface in ABAQUS CAE

3. plugin3_plugin.py -
This code is responsible for creating the plugin button and its logo within the ABAQUS CAE menu

An instructional video for the use of this tool is found at https://youtu.be/28rQcXDkIVA
