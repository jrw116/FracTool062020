Plug-in tool within ABAQUS CAE to perform energy and stength calculations along any given desired fracture plane

There are three main files within the plugin:

1. EDC2.py -
This is the main file responsible for performing all calculations to determine the energy and strength of desired tensile fracture planes.
	- Within this file there are multiple functions that are used:
		- Plycomp: if a hexahedral mesh is chosen this finds the cut plane area in each element split by the number of plies
		- Plycomptet: if a tetrahedral mesh is chosen this finds the cut plane area in each element split by the number of plies
		- Plycomphextet: if a mixed hexahedral and tetrahedral mesh is chosen this finds the cut plane area in each element split by the number of plies
		- Plycomphexel: for a mixed mesh this find the cut plane area in each element for hexahedral elements
		- Plycomptetel: for a mixed mesh this find the cut plane area in each element for tetrahedral elements 
		- Line: sorts points for the tetrahedral mesh
		- Poly4: calculates area of quadrilateral with 3D nodes
		- Poly5: calculates area of polygons with more than 4 3D nodes
		- FibreMat: calculates failure energy and force for given plane
		- Plot: plots graphs and creates Results.csv file with results of running the plug-in

2. plugin3DB.py - 
This file is responsible for creating the user interface in ABAQUS CAE

3. plugin3_plugin.py -
This file is responsible for creating the plugin button and its logo within the ABAQUS CAE menu

An instructional video for the use of this tool is found at https://youtu.be/28rQcXDkIVA
