# minedgefactory
Source code for a mixed-integer linear programming (MILP) approach to solve the minimum-hyperedge factory problem.
The main file is run.py, which will find the minimum-hyperedge factory under conservation and accumulation for no negative regulation and 1st order negative regulation, is run with the following usage:

For sbml files (which must have already been parsed by xmlparser.py):
`python run.py --name <name of -hypernodes.txt file> --sbml --target <desired target>`

For Reactome with stoichiometry:
`python run.py --name <name of -hypernodes.txt file> --stoichiometry --target <desired target>`

For second order negative regulation on Reactome:
`python run.py --name <name of -hypernodes.txt file> --secondordernegreg --target <desired target>`

Under construction...
