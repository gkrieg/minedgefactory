import xml.etree.ElementTree as ET
import sys
tree = ET.parse(sys.argv[1])
root = tree.getroot()
#print(root[0][6][0].attrib)
rootdir = '../biopax-hypergraph/test/parsed/'
#rootdir = ''
outf = open('{}parsed{}-hyperedges.txt'.format(rootdir,sys.argv[2]),'w')
nodesoutf = open('{}parsed{}-hypernodes.txt'.format(rootdir,sys.argv[2]),'w')

nodes = set()
outf.write('#Tail Head PosReg NegReg ID Stoichiometries\n')
for reaction in root[0][6]:
    isreversible = True if 'true' in reaction.attrib['reversible'] else False
    for field in reaction:
        if 'listOfReactants' in field.tag:
            reactantnodes = field 
        elif 'listOfProducts' in field.tag:
            productnodes = field 
    reactants = [r.attrib['species'] for r in reactantnodes]
    rstoi = ['{}:{}'.format(r.attrib['species'],r.attrib['stoichiometry']) for r in reactantnodes]
    products = [r.attrib['species'] for r in productnodes]
    pstoi = ['{}:{}'.format(r.attrib['species'],r.attrib['stoichiometry']) for r in productnodes]
    for r in reactants:
        nodes.add(r)
    for r in products:
        nodes.add(r)
    outf.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(';'.join(reactants),';'.join(products),'None', 'None', 'Nullstring',';'.join(rstoi + pstoi)))
    if isreversible == True:
        outf.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(';'.join(products),';'.join(reactants),'None', 'None', 'Nullstring',';'.join(rstoi + pstoi)))

nodesoutf.write('#hypernode nodes\n')
for n in nodes:
    nodesoutf.write('{}\t{}\n'.format(n,n))


