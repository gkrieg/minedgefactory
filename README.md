# minedgefactory
Source code for a mixed-integer linear programming (MILP) approach to solve the minimum-hyperedge factory problem.

The python code requires installation of the Halp package in Python 2.7, and CPLEX.
In the datasets folder, the two hypernodes.txt files for allreactomestoichiometry must be concatenated together before running that dataset.

The main file is run.py, which will find the minimum-hyperedge factory under conservation and accumulation for no negative regulation and 1st order negative regulation, is run with the following usage:

For sbml files (which must have already been parsed by xmlparser.py):
`python run.py --name <name of -hypernodes.txt file> --sbml --target <desired target>`

For Reactome with stoichiometry:
`python run.py --name <name of -hypernodes.txt file> --stoichiometry --target <desired target>`

For first order negative regulation on Reactome:
`python run.py --name <name of -hypernodes.txt file> --stoichiometry --negative --target <desired target>`

For second order negative regulation on Reactome:
`python run.py --name <name of -hypernodes.txt file> --stoichiometry --secondordernegreg --target <desired target>`

The repository also comes with a parser that takes SBML files and converts them into the hyperedges.txt and hypernodes.txt files required for the application.
To run the parser:
`python xmlparser.py <sbml file name>`

Under construction...
