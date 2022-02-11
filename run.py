#!/usr/bin/python

## ISMB poster submission
import os
import os.path
import sys
from optparse import OptionParser
import time
import pickle as pkl
import traceback

# for plotting

from stoichiometryilp import *

from halp.directed_hypergraph import DirectedHypergraph
from halp.utilities import directed_graph_transformations
from halp.algorithms.directed_paths import b_visit

## http://manual.graphspace.org/en/latest/Programmers_Guide.html#graphspace-python
## http://manual.graphspace.org/projects/graphspace-python/en/latest/
#from graphspace_python.graphs.classes.gsgraph import GSGraph
#from graphspace_python.api.client import GraphSpace
#import graphspace_interface as interface


## GLOBAL VARIABLES
DELIM=';'
FORCE = False
PRINTONLY = False
ROOTDIR = '/Users/skrieger/Documents/UofA/biopax-hypergraph/biopax-hypergraph/test'

COLORS = {'blue':'#6F95EB',
        'red':'#EB856F',
        'orange':'#EBCD6F',
        'gray':'#C6C6C6',
        'black':'#000000',
        'white':'#FFFFFF'
        }

def main(args):
    opts = parseOptions(args)

    if not opts.stoichiometry:
        H,_,__ = make_hypergraph(ROOTDIR+'/parsed/{0}'.format(opts.name))
        print(('target: {0}'.format(opts.target)))
        H_sources,H_targets,H_high_penalty_sources,H_name_dict = getSourcesTargets2(opts.name,H,'hypergraph',opts.source,opts.target)
        source,target = add_super_nodes(H,H_sources,H_targets,H_high_penalty_sources,opts.name)

    # run
    if opts.minimal_precursor:
        H,_,__,stoichiometries,negative_regulators = make_hypergraph(ROOTDIR+'/parsed/{0}'.format(opts.name),stoichiometry=True,keep_negative_regulators=True)
        #print('numhyperedges',len(H.get_hyperedge_id_set()))
        #print(len(stoichiometries))
        #print(stoichiometries,negative_regulators)
        print(('target: {0}'.format(opts.target)))
        sources,targets,high_penalty_sources,H_name_dict = getSourcesTargets2(opts.name,H,'hypergraph',opts.source,opts.target)
        sources.update(high_penalty_sources)
        print('len sources',len(sources))
        #source,target,stoichiometries = add_super_nodes(H,H_sources,H_targets,H_high_penalty_sources,opts.name,stoichiometries=stoichiometries)
        ilpname = '%s/results/reactome-%s-%s-precursors.lp' % (ROOTDIR,opts.name,opts.type)
        outprefix = '%s/results/reactome-%s-%s-precursors' % (ROOTDIR,opts.name,opts.type)
        #source_edges,H,stoichiometries,negative_regulators = find_source_edges(H,sources,add_edges=True,stoichiometries=stoichiometries,negative_regulators=negative_regulators)
        source_edges = find_source_edges(H,sources,add_edges=False,stoichiometries=stoichiometries,negative_regulators=negative_regulators)
        ilpfile = open(ilpname,'w')
        write_binary = True

        write_minimal_precursor_ilp(ilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary)
        ilpfile.close()
        try:
            solve_minimal_precursor_ilp(ilpname,H,binary_vars=write_binary)
            print('solved first precursor MILP\n')
            #print('neg\n')
        except:
            print('not able to solve first precursor ILP\n')
        ilpfile = open(ilpname,'w')
        write_binary = True
        write_minimal_precursor_ilp(ilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=False,binary_vars=write_binary)
        ilpfile.close()
        try:
            solve_minimal_precursor_ilp(ilpname,H,binary_vars=write_binary)
            print('solved precursor conservation MILP\n')
        except:
            print('not able to solve precursor conservation ILP\n')


        ilpfile = open(ilpname,'w')
        write_binary = True
        write_minimal_precursor_ilp(ilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary,negative_regulators=negative_regulators)
        ilpfile.close()
        try:
            solve_minimal_precursor_ilp(ilpname,H,binary_vars=write_binary)
            print('solved first negative precursor MILP\n')
            #print('neg\n')
        except:
            print('not able to solve first negative precursor ILP\n')
        ilpfile = open(ilpname,'w')
        write_binary = True
        write_minimal_precursor_ilp(ilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=False,binary_vars=write_binary,negative_regulators=negative_regulators)
        ilpfile.close()
        try:
            solve_minimal_precursor_ilp(ilpname,H,binary_vars=write_binary)
            print('solved negative precursor conservation MILP\n')
        except:
            print('not able to solve negative precursor conservation ILP\n')

    if opts.hybrid:
        H,_,__,stoichiometries,reversible_reactions = make_hypergraph(ROOTDIR+'/parsed/{0}'.format(opts.name),stoichiometry=True,return_reversible=True)
        print(('target: {0}'.format(opts.target)))
        sources,targets,high_penalty_sources,H_name_dict = getSourcesTargets2(opts.name,H,'hypergraph',opts.source,opts.target)
        sources.update(high_penalty_sources)
        print('len sources',len(sources))
        pilpname = '%s/results/reactome-%s-%s-precursors.lp' % (ROOTDIR,opts.name,opts.type)
        poutprefix = '%s/results/reactome-%s-%s-precursors' % (ROOTDIR,opts.name,opts.type)

        milpname = '%s/results/reactome-%s-%s-min.lp' % (ROOTDIR,opts.name,opts.type)
        moutprefix = '%s/results/reactome-%s-%s-min' % (ROOTDIR,opts.name,opts.type)

        silpname = '%s/results/reactome-%s-%s-acyclic.lp' % (ROOTDIR,opts.name,opts.type)
        soutprefix = '%s/results/reactome-%s-%s-acyclic' % (ROOTDIR,opts.name,opts.type)
        #source_edges,H,stoichiometries,negative_regulators = find_source_edges(H,sources,add_edges=True,stoichiometries=stoichiometries,negative_regulators=negative_regulators)
        source_edges = find_source_edges(H,sources,add_edges=False,stoichiometries=stoichiometries)
        write_binary = True
        begintime = time.time()
        ones = []

        milpfile = open(milpname,'w')

        write_stoichiometry_ilp(milpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary,sourcevars=True)
        milpfile.close()
        try:
            __, ones = solve_stoichiometry_ilp(milpname,H,sources,binary_vars=write_binary)
            print('solved first MILP\n')
            #print('neg\n')
        except:
            print('not able to solve first ILP\n')
        sftime1 = time.time()
        print('time for first sf MILP {}'.format(sftime1-begintime))

        if ones != []:

            pilpfile = open(pilpname,'w')
            write_minimal_precursor_ilp(pilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary,reversible=reversible_reactions,edgenum=len(ones))
            #write_minimal_precursor_ilp(pilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary,reversible=False)
            pilpfile.close()
            sourceones = []
            try:
                _, sourceones = solve_minimal_precursor_ilp(pilpname,H,binary_vars=write_binary)
                print('solved first precursor MILP\n')
                #print('neg\n')
            except:
                print('not able to solve first precursor ILP\n')
            minpretime1 = time.time()
            print('time for first mps MILP {}'.format(minpretime1-sftime1))
        else:
            print('did not attempt to solve mpsMILP')

            

    if opts.sbml:
        H,_,__,stoichiometries,reversible_reactions = make_hypergraph(ROOTDIR+'/parsed/{0}'.format(opts.name),stoichiometry=True,return_reversible=True)
        print(('target: {0}'.format(opts.target)))
        sources,targets,high_penalty_sources,H_name_dict = getSourcesTargets2(opts.name,H,'hypergraph',opts.source,opts.target)
        sources.update(high_penalty_sources)
        print('len sources',len(sources))
        pilpname = '%s/results/reactome-%s-%s-precursors.lp' % (ROOTDIR,opts.name,opts.type)
        poutprefix = '%s/results/reactome-%s-%s-precursors' % (ROOTDIR,opts.name,opts.type)

        milpname = '%s/results/reactome-%s-%s-min.lp' % (ROOTDIR,opts.name,opts.type)
        moutprefix = '%s/results/reactome-%s-%s-min' % (ROOTDIR,opts.name,opts.type)

        silpname = '%s/results/reactome-%s-%s-acyclic.lp' % (ROOTDIR,opts.name,opts.type)
        soutprefix = '%s/results/reactome-%s-%s-acyclic' % (ROOTDIR,opts.name,opts.type)
        #source_edges,H,stoichiometries,negative_regulators = find_source_edges(H,sources,add_edges=True,stoichiometries=stoichiometries,negative_regulators=negative_regulators)
        source_edges = find_source_edges(H,sources,add_edges=False,stoichiometries=stoichiometries)
        write_binary = True
        begintime = time.time()

        pilpfile = open(pilpname,'w')
        write_minimal_precursor_ilp(pilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary,reversible=reversible_reactions)
        #write_minimal_precursor_ilp(pilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary,reversible=False)
        pilpfile.close()
        try:
            solve_minimal_precursor_ilp(pilpname,H,binary_vars=write_binary)
            print('solved first precursor MILP\n')
            #print('neg\n')
        except:
            print('not able to solve first precursor ILP\n')
        minpretime1 = time.time()
        print('time for first precursor MILP {}'.format(minpretime1-begintime))
        pilpfile = open(pilpname,'w')
        write_minimal_precursor_ilp(pilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=False,binary_vars=write_binary,reversible=reversible_reactions)
        #write_minimal_precursor_ilp(pilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=False,binary_vars=write_binary,reversible=False)
        pilpfile.close()
        try:
            solve_minimal_precursor_ilp(pilpname,H,binary_vars=write_binary)
            print('solved precursor conservation MILP\n')
        except:
            print('not able to solve precursor conservation ILP\n')
        minpretime2 = time.time()
        print('time for conservation precursor MILP {}'.format(minpretime2-minpretime1))
        milpfile = open(milpname,'w')

        write_stoichiometry_ilp(milpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary,reversible=reversible_reactions)
        #write_stoichiometry_ilp(milpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary,reversible=False)
        milpfile.close()
        try:
            solve_stoichiometry_ilp(milpname,H,sources,binary_vars=write_binary)
            print('solved first MILP\n')
            #print('neg\n')
        except:
            print('not able to solve first ILP\n')
        sftime1 = time.time()
        print('time for first sf MILP {}'.format(sftime1-minpretime2))

        milpfile = open(milpname,'w')
        write_stoichiometry_ilp(milpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=False,binary_vars=write_binary,reversible=reversible_reactions)
        #write_stoichiometry_ilp(milpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=False,binary_vars=write_binary,reversible=False)
        milpfile.close()
        try:
            solve_stoichiometry_ilp(milpname,H,sources,binary_vars=write_binary)
            print('solved conservation MILP\n')
        except:
            print('not able to solve conservation ILP\n')
        sftime2 = time.time()
        print('time for conservation sf MILP {}'.format(sftime2-sftime1))



        source,target = add_super_nodes(H,sources,targets,high_penalty_sources,opts.name)
        H, vertex_dict = convert_hypergraph_nodes(H)
        ilpname = '%s/results/reactome-%s-%s-acyclic.lp' % (ROOTDIR,opts.name,opts.type)
        outprefix = '%s/results/reactome-%s-%s-acyclic' % (ROOTDIR,opts.name,opts.type)
        make_shortest_acyclic_hyperpath_ilp(H,source,target,silpname)


        variables,objective = runILP_singlesol(H,silpname,soutprefix,opts.force,target,source)
        print('True objective is ', objective-1)
        if objective > -1:
            pathedges = [v for v in variables[0] if variables[0][v] > 0.2 and 'e' in v]
            for v in variables[0]:
                if variables[0][v] > 0.2 and 'e' in v:
                    print(v,variables[0][v])
                    #print(H.get_hyperedge_tail(v[2:]),H.get_hyperedge_head(v[2:]))
        sptime = time.time()
        print('time for shortest path MILP {}'.format(sptime-sftime2))

    if opts.stoichiometry:
        timebegin = time.time()
        H,_,__,stoichiometries,negative_regulators = make_hypergraph(ROOTDIR+'/parsed/{0}'.format(opts.name),stoichiometry=True,keep_negative_regulators=True)
        #print('numhyperedges',len(H.get_hyperedge_id_set()))
        #print(len(stoichiometries))
        #print(stoichiometries,negative_regulators)
        print(('target: {0}'.format(opts.target)))
        sources,targets,high_penalty_sources,H_name_dict = getSourcesTargets2(opts.name,H,'hypergraph',opts.source,opts.target)
        sources.update(high_penalty_sources)
        print('len sources',len(sources))
        #source,target,stoichiometries = add_super_nodes(H,H_sources,H_targets,H_high_penalty_sources,opts.name,stoichiometries=stoichiometries)
        ilpname = '%s/results/reactome-%s-%s-stoichiometry2.lp' % (ROOTDIR,opts.name,opts.type)
        outprefix = '%s/results/reactome-%s-%s-stoichiometry2' % (ROOTDIR,opts.name,opts.type)
        #source_edges,H,stoichiometries,negative_regulators = find_source_edges(H,sources,add_edges=True,stoichiometries=stoichiometries,negative_regulators=negative_regulators)
        source_edges = find_source_edges(H,sources,add_edges=False,stoichiometries=stoichiometries,negative_regulators=negative_regulators)
        if opts.second_order_neg_reg:
            timebegin = time.time()
            solution,sourcesused,i = second_order_negative_regulation(ilpname,H,stoichiometries,targets,source_edges,negative_regulators,sources=sources,positiveflux=True)
            print(solution,len(solution))
            print(sourcesused,len(sourcesused))
            print('{} iterations needed'.format(i))
            if solution != '':
                print('solved second order neg reg')
            else:
                print('failed to find second order neg reg solution')

            timeaccend = time.time()
            print('time for second order negative regulation {}'.format(timeaccend-timebegin))

            solution,sourcesused,i = second_order_negative_regulation(ilpname,H,stoichiometries,targets,source_edges,negative_regulators,sources=sources,positiveflux=False)
            print(solution,len(solution))
            print(sourcesused,len(sourcesused))
            print('{} iterations needed'.format(i))
            if solution != '':
                print('solved second order neg reg conservation')
            else:
                print('failed to find second order neg reg solution conservation')

            timeconsend = time.time()
            print('time for second order negative regulation {}'.format(timeconsend-timeaccend))

        elif opts.facenumeration:

            timebegin = time.time()
            solution,sourcesused,i = enumerate_minimal_factories(ilpname,H,stoichiometries,targets,source_edges,negative_regulators,sources=sources,positiveflux=True)
            print(solution,len(solution))
            print(sourcesused,len(sourcesused))
            print('{} iterations needed'.format(i))
            timeend = time.time()
            print('time taken is {}'.format(timeend-timebegin))

        else:
            timebegin = time.time()
            ilpfile = open(ilpname,'w')
            write_binary = True
            write_stoichiometry_ilp(ilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary)
            #write_stoichiometry_ilp(ilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary,negative_regulators=negative_regulators)
            ilpfile.close()
            try:
                solve_stoichiometry_ilp(ilpname,H,sources,binary_vars=write_binary,solnum=3)
                print('solved first MILP\n')
                #print('neg\n')
            except:
                
                print('not able to solve first ILP\n')
            ilpfile = open(ilpname,'w')
            write_binary = True
            write_stoichiometry_ilp(ilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=False,binary_vars=write_binary)
            #write_stoichiometry_ilp(ilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=False,binary_vars=write_binary,negative_regulators=negative_regulators)
            ilpfile.close()
            try:
                solve_stoichiometry_ilp(ilpname,H,sources,binary_vars=write_binary,solnum=3)
                print('solved conservation MILP\n')
            except:
                print('not able to solve conservation ILP\n')

            if opts.negative:

                ilpfile = open(ilpname,'w')
                write_binary = True
                #write_stoichiometry_ilp(ilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary)
                write_stoichiometry_ilp(ilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=True,binary_vars=write_binary,negative_regulators=negative_regulators)
                ilpfile.close()
                try:
                    solve_stoichiometry_ilp(ilpname,H,sources,binary_vars=write_binary,solnum=3)
                    print('solved first MILP negative regulation\n')
                    #print('neg\n')
                except:
                    
                    print('not able to solve first ILP negative regulation\n')
                ilpfile = open(ilpname,'w')
                write_binary = True
                #write_stoichiometry_ilp(ilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=False,binary_vars=write_binary)
                write_stoichiometry_ilp(ilpfile,H,stoichiometries,targets,source_edges,sources=sources,positiveflux=False,binary_vars=write_binary,negative_regulators=negative_regulators)
                ilpfile.close()
                try:
                    solve_stoichiometry_ilp(ilpname,H,sources,binary_vars=write_binary,solnum=3)
                    print('solved conservation MILP negative regulation\n')
                except:
                    print('not able to solve conservation ILP negative regulation\n')

            #source,target = add_super_nodes(H,sources,targets,high_penalty_sources,opts.name)
            #H, vertex_dict = convert_hypergraph_nodes(H)
            #done added by Spencer
            #ilpname = '%s/results/reactome-%s-%s-acyclic.lp' % (ROOTDIR,opts.name,opts.type)
            #outprefix = '%s/results/reactome-%s-%s-acyclic' % (ROOTDIR,opts.name,opts.type)
            #make_shortest_acyclic_hyperpath_ilp(H,source,target,ilpname)


            #where time bound code goes
            #variables,objective = runILP_singlesol(H,ilpname,outprefix,opts.force,target,source)
            #print('True objective is ', objective-1)
            #if objective > -1:
                #pathedges = [v for v in variables[0] if variables[0][v] > 0.2 and 'e' in v]
                #for v in variables[0]:
                    #if variables[0][v] > 0.2 and 'e' in v:
                        #print(v,variables[0][v])
                        #print(H.get_hyperedge_tail(v[2:]),H.get_hyperedge_head(v[2:]))
            
            #node_dict = pkl.load(open('node_dict.pkl','rb'))
            #taildistancelist,H2,P,tailpathlist = tail_path_heuristic(H,source,target,node_dict=node_dict,return_paths=True)
            #print(len(P),P)
            #make_shortest_hyperpath_ilp_simple(H,source,target,ilpname)
            #variables = runILP_singlesol(H,ilpname,outprefix,opts.force,target)
            


def convert_hypergraph_nodes(H):
    H2 = DirectedHypergraph()
    edgecount = len(H.get_hyperedge_id_set())
    edgenums = sorted(int(e[1:]) for e in H.get_hyperedge_id_set())
    vertexdict = {}
    nodenum = 1
    for node in H.node_iterator():
        if node != 'SUPERTARGET' and node != 'SUPERSOURCE':
            vertexdict[node] = nodenum
            nodenum += 1
        else:
            vertexdict[node] = node
    for e in edgenums:
        if H.has_hyperedge_id('e{0}'.format(e)):
            tail = [vertexdict[v] for v in H.get_hyperedge_tail('e{0}'.format(e))]
            for node in tail:
                if type(node) is not int:
                    print(node)
            head = [vertexdict[v] for v in H.get_hyperedge_head('e{0}'.format(e))]
            for node in head:
                if type(node) is not int:
                    print(node)
            H2.add_hyperedge(set(tail),set(head),weight=H.get_hyperedge_weight('e{0}'.format(e)))
    for node in H2.node_iterator():
        if type(node) is not int:
            print(node)
    for edge in H2.hyperedge_id_iterator():
        for node in H2.get_hyperedge_tail(edge):
            if type(node) is not int:
                print(node)
        for node in H2.get_hyperedge_head(edge):
            if type(node) is not int:
                print(node)
    return H2,vertexdict

def runILP_singlesol(H,ilpname,outprefix,force,target,source,cheats=None):
    objective=-1
    if force or not os.path.isfile('%s-1.variables' % (outprefix)):
        if not cheats:
            numsols,numoptobjective,allvars,objective = solveILP(H,H.get_node_set(),ilpname,outprefix,1,target,source,printobjective=True)
        else:
            numsols,numoptobjective,allvars = solveILP_cheats(H,H.get_node_set(),ilpname,outprefix,1,cheats)
        #allvars = allvars[0]
    else:
        print('not solving ILP. Use --force to override.')
        allvars = {}
        with open('%s-1.variables' % (outprefix)) as fin:
            for line in fin:
                if line[0] == '#': 
                    continue
                row = line.strip().split()
                if 'o_' in row[0]:
                    allvars[row[0]] = float(row[1])
                else:
                    allvars[row[0]] = int(row[2])
    return allvars,objective
#############################################

def parseOptions(args):
    desc = 'python master-script.py [options]'
    parser = OptionParser(usage=desc)

    # General Options
    parser.add_option('','--force',action='store_true',help='Overwrite files if they exist.')
    parser.add_option('','--printonly',action='store_true',help='Print commands to screen, but do not execute them.')
    
    # EXPERIMENTS/TESTS
    parser.add_option('','--name',type='string',default='WNT5A',help='Name of dataset (WNT5A, CTNNB1, WNT, or ALL). Default=WNT.')
    parser.add_option('','--type',type='string',default='hypergraph',help='graph type: hypergraph, graph-with-complexes, or graph. Default=hypergraph.')
    parser.add_option('','--stoichiometry',action='store_true',help='stoichiometric hyperpath from S to T.')
    parser.add_option('','--minimal_precursor',action='store_true',help='stoichiometric hyperpath from S to T.')
    parser.add_option('','--sbml',action='store_true',help='stoichiometric hyperpath from S to T.')
    parser.add_option('','--negative',action='store_true',help='stoichiometric hyperpath from S to T with first order negative regulation.')
    parser.add_option('','--second_order_neg_reg',action='store_true',help='stoichiometric hyperpath from S to T.')
    parser.add_option('','--hybrid',action='store_true',help='Compute heuristic shortest cyclic hyperpath from S to T.')
    parser.add_option('','--facenumeration',action='store_true',help='enumerate heuristic shortest cyclic hyperpath from S to T.')
    parser.add_option('','--source',type='string',action='append',help='Sources. Default = WNT5A')
    parser.add_option('','--target',type='string',action='append',help='Targets. Default = CTNNB1')

    # VISUALIZE
    parser.add_option('','--viz',action='store_true',help='Post to GraphSpace')

    opts,args = parser.parse_args()
       ## set FORCE and PRINTONLY variables
    global FORCE, PRINTONLY
    if opts.force:
       FORCE = True
    if opts.printonly:
       PRINTONLY = True

    if not opts.source:
       opts.source = ['WNT5A']
    if not opts.target:
       opts.target = ['CTNNB1']
    print('OPTIONS ARE',opts)
    return opts


def make_hypergraph(file_prefix,delim=';',sep='\t',keep_singleton_nodes=False,stoichiometry=False,keep_negative_regulators=False,return_reversible=False):
    hypernodes = {}
    with open(file_prefix+'-hypernodes.txt') as fin:
        for line in fin:
            if line[0] == '#':
                continue
            row = line.strip().split(sep)
            ## TODO fix -- this happens when a complex contains other complexes
            if len(row) == 1:
                hypernodes[row[0]] = ['OtherComplexes-FIX']
            else:
                hypernodes[row[0]] = row[1].split(delim)
    print(('%d hypernodes from hypernodes file' % (len(hypernodes))))
    identifier2id = {}
    id2identifier = {}
    H = DirectedHypergraph()
    if keep_singleton_nodes:
        for n in hypernodes:
            H.add_node(n)

    skipped1 = 0
    skipped2 = 0
    tailsizes = []
    headsizes = []
    selfloops = []
    noinselfloops = 0
    indegree = []
    outdegree = []
    numtargets = 0
    numsources = 0
    if stoichiometry == True:
        stoichiometries = {}
    if keep_negative_regulators == True:
        negative_regulators = {}
    if return_reversible == True:
        reversible_reactions = {}
        prevhead = set()
        prevtail = set()
        previd = None

    with open(file_prefix+'-hyperedges.txt') as fin:
        for line in fin:
            if line[0] == '#':
                continue
            #row = line.strip().split(sep)
            row = line.strip().split()
            tail = set()
            head = set()

            ## Tail includes tail and regulators.
            ## Head includes head.
            if row[0] != 'None' and row[0] != '':
                tail.update(row[0].split(delim))
            if row[1] != 'None' and row[1] != '':
                head.update(row[1].split(delim))
            if row[2] != 'None' and row[2] != '':
                tail.update(row[2].split(delim))
            #These are the negative regulators!
            #if row[3] != 'None':
                #tail.update(row[3].split(delim))
            hedge_id = row[4]

            ## THIS IS A HACK FOR NOW ( should be incorporated in the make-hypergraph.py code)
            ## IGnore any reactions that have a Reactome Identifier (e.g. has "HSA") instead of
            ## a PAthway Commons identifier.
            if any(['HSA' in s for s in tail]+['HSA' in s for s in head]):
                skipped1+=1
            elif len(tail)==0 or len(head)==0:
                skipped2+=1
            else:
                hid = H.add_hyperedge(tail,head,identifier=hedge_id)
                tailsizes.append(len(tail))
                headsizes.append(len(head))
                intersection = tail.intersection(head)
                if len(intersection) > 0:
                    selfloops.append([v for v in intersection])

                identifier2id[hedge_id] = hid
                id2identifier[hid] = hedge_id
                if keep_negative_regulators == True:
                    if row[3] != 'None':
                        negative_regulators[hid] = row[3].split(delim)
                if stoichiometry == True:
                    stoichiometries[hid] = {}
                    if row[5] != 'None':
                        for molstoi in row[5].split(delim):
                            molecule = 'http:{}'.format(molstoi.split(':')[1])
                            if 'parsed' in file_prefix.split('/')[-1]:
                                stoichio = molstoi.split(':')[1]
                            else:
                                stoichio = molstoi.split(':')[2]
                            stoichiometries[hid][molecule] = float(stoichio)
                    for molecule in tail:
                        if molecule not in stoichiometries[hid]:
                            stoichiometries[hid][molecule] = 1.0
                    for molecule in head:
                        if molecule not in stoichiometries[hid]:
                            stoichiometries[hid][molecule] = 1.0
                if return_reversible == True:
                    if head == prevtail and tail == prevhead:
                        reversible_reactions[hid] = previd
                        reversible_reactions[previd] = hid
                    prevhead = head
                    prevtail = tail
                    previd = hid

    if 'WNT' in file_prefix.split('/') :
        print(('num vertices {0}'.format(len(H.get_node_set()))))
        print(('num edges {0}'.format(len(H.get_hyperedge_id_set()))))
        print(('largest head {0}'.format(max(headsizes))))
        print(('median head {0}'.format(median(headsizes))))
        print(('largest tail {0}'.format(max(tailsizes))))
        print(('median tail {0}'.format(median(tailsizes))))
        print(('num selfloops {0}'.format(len(selfloops))))
        selfset = set()
        for l in selfloops:
            selfset.update(l)
        for v in selfset:
            intersection = H.get_backward_star(v).intersection(H.get_forward_star(v))
            if len(intersection) == len(H.get_backward_star(v)):
                noinselfloops += 1
        print(('num no in selfloops {0}'.format(noinselfloops)))
        for v in H.get_node_set():
            outdegree.append(len(H.get_forward_star(v)))
            indegree.append(len(H.get_backward_star(v)))
            if len(H.get_forward_star(v)) == 0:
                numtargets += 1
            if len(H.get_backward_star(v)) == 0:
                numsources += 1
        print(('max outdegree {0}'.format(max(outdegree))))
        print(('max indegree {0}'.format(max(indegree))))
        print(('median indegree {0}'.format(median(indegree))))
        print(('median outdegree {0}'.format(median(outdegree))))
        print(('num targets {0}'.format(numtargets)))
        print(('num sources {0}'.format(numsources)))
    print(('%d reactions skipped because of Reactome identifier' % (skipped1)))
    print(('%d reactions skipped because of an empty tail or head' % (skipped2)))
    ## annotate nodes
    num_hypernodes = 0
    for node in H.get_node_set():
        if node in hypernodes and hypernodes[node] != [node]:
            H.add_node(node,hypernode_members=hypernodes[node],is_hypernode=True)
            num_hypernodes+=1
        else:
            H.add_node(node,is_hypernode=False,hypernode_members=[])

        H.add_node(node)

    #print('Hypergraph has %d hyperedges and %d nodes (%d of these are hypernodes)' %
        #(stats.number_of_hyperedges(H),stats.number_of_nodes(H),num_hypernodes))

    if stoichiometry == False:
        if keep_negative_regulators == False:
            if return_reversible == False:
                return H, identifier2id, id2identifier
            else:
                return H, identifier2id, id2identifier, reversible_reactions

        else:
            if return_reversible == False:
                return H, identifier2id, id2identifier, negative_regulators
            else:
                return H, identifier2id, id2identifier, negative_regulators, reversible_reactions
    else:
        if keep_negative_regulators == False:
            if return_reversible == False:
                return H, identifier2id, id2identifier, stoichiometries
            else:
                return H, identifier2id, id2identifier, stoichiometries, reversible_reactions
        else:
            if return_reversible == False:
                return H, identifier2id, id2identifier, stoichiometries, negative_regulators
            else:
                return H, identifier2id, id2identifier, stoichiometries, negative_regulators, reversible_reactions
        

def getSourcesTargets2(name,H,graph_type,source_list,target_list):
    if name == 'WNT5A':
        sources = set(['http://pathwaycommons.org/pc12/Complex_850f44f917acb11059ff27c88a0494ee'])
        #sources = set(['http://pathwaycommons.org/pc12/Protein_037001d14ad7601b82c325eeaac1cc36']) #example 2 to break the loop
        #sources = set(['http://pathwaycommons.org/pc12/Protein_bdf326eedb65f7fe6a0259c5cb8c4ed4','http://pathwaycommons.org/pc12/Protein_037001d14ad7601b82c325eeaac1cc36']) #Trying to find cyclic

        targets = set(['http://pathwaycommons.org/pc12/Protein_6eae9e1fb8e906b20a3ebdf4485a4a3d']) #example 1
        #targets = set(['http://pathwaycommons.org/pc12/Protein_bc47c96d7c652d22f94260b30d5c8043']) #example 2
        #targets = set(['http://pathwaycommons.org/pc12/Complex_46f99a13ca1b39c8d93926b9b394c395'])# trying to find cyclic
    elif name == 'WNT':
        sources = set(['http://pathwaycommons.org/pc12/Complex_850f44f917acb11059ff27c88a0494ee','http://pathwaycommons.org/pc12/Protein_355a3029f445775a6c82451d5c86031b','http://pathwaycommons.org/pc12/Protein_6b903a8964fcd5874a866c1e38e37456'])
        targets = set(['http://pathwaycommons.org/pc12/Protein_3c9a7ce5eec5c6ec97c1a3008c3c9c99','http://pathwaycommons.org/pc12/Protein_1625818ba862a465e6bfe45c1a57c0ec'])
    elif name == 'allpid':
        sources = set(['http://pathwaycommons.org/pc12/Complex_850f44f917acb11059ff27c88a0494ee','http://pathwaycommons.org/pc12/Protein_355a3029f445775a6c82451d5c86031b','http://pathwaycommons.org/pc12/Protein_6b903a8964fcd5874a866c1e38e37456'])
        targets = set(['http://pathwaycommons.org/pc12/Complex_81ba3b0707b6c6abd477dd59315147f4'])
        if ':' in target_list[0]:
            targets = set([target_list[0]])
    elif name == 'allreactome' or name == 'allreactomestoichiometry':
        sources = set(['http://pathwaycommons.org/pc12/Complex_eb845eb263044d3b435b479eb76ac674'])
        if len(target_list) < 1:
        #targets = set(['http://pathwaycommons.org/pc12/Protein_83baeb7dc5ecdcda2877b086aebb603f'])
            targets = set(['http://pathwaycommons.org/pc12/Complex_2f0a923dd4cf0b828c8176a962e14011'])
        else:
            targets = set(target_list)
    elif name == 'WNTSTOI2':
        sources = set(['s'])
        targets = set(['t'])
    elif name == 'WNTSTOI':
        sources = set()
        targets = set([v for v in H.node_iterator() if len(H.get_forward_star(v)) == 0])
    elif len(target_list) > 0:
        sources = set()
        targets = set(target_list)
    else:
        sources = set()
        targets = set([min(H.get_node_set())])

    print(('sources,targets',sources,targets))
    for v in targets:
        for e in H.get_backward_star(v):
            print((e,H.get_hyperedge_tail(e),H.get_hyperedge_head(e)))
    #for e in H.hyperedge_id_iterator():
        #for v in H.get_hyperedge_head(e):
            #if v in targets:
                #print(e,H.get_hyperedge_tail(e),H.get_hyperedge_head(e))

    ## backward star sources
    high_penalty_sources = set([n for n in H.get_node_set() if len(H.get_backward_star(n))==0]).difference(sources)
    #all_sources = set(sources.union(high_penalty_sources))
    all_sources = sources 
    if len(all_sources.intersection(targets)) >0:
        print('Warning: removing %d nodes from targets that were both sources and targets' % (len(set(sources).intersection(set(targets)))))
        targets = [t for t in targets if t not in all_sources]

    print('%d sources and %d targets; %d high penalty sources (empty backward stars)'  % (len(sources),len(targets),len(high_penalty_sources)))
    name_dict = {}
    return sources,targets,high_penalty_sources,name_dict


def add_super_nodes(H,sources,targets,high_penalty_sources,dataset,stoichiometries=None):

    super_source = 'SUPERSOURCE'
    print('-------------------------------sourcelen is: ',len(sources))
    largenum = len(H.get_hyperedge_id_set())+100
    #for s in sources:
        #H.add_hyperedge(set([super_source]), set([s]), weight=1)
    hid = H.add_hyperedge(set([super_source]),set(sources.union(high_penalty_sources)),weight=0)
    if stoichiometries != None:
        stoichiometries[hid] = {}
        stoichiometries[hid][super_source] = 1.0
        for s in set(sources.union(high_penalty_sources)):
            stoichiometries[hid][s] = 1.0
    
    #for s in high_penalty_sources:
    ##this is where I changed the weight of high-penalty sources -Anna
        #H.add_hyperedge(set([super_source]), set([s]), weight=1)    
        #H.add_hyperedge(set([super_source]), set([s]), weight=100)    

    super_target = 'SUPERTARGET'
    if dataset == 'test' or dataset == 'reactome' or dataset == 'WNT':
        hid = H.add_hyperedge(set(targets),set([super_target]), weight=1)
        if stoichiometries != None:
            stoichiometries[hid] = {}
            stoichiometries[hid][super_target] = 1.0
            for s in set(targets):
                stoichiometries[hid][s] = 1.0
    else: # dataset == 'ncipid'
        for t in targets:
            hid = H.add_hyperedge(set([t]),set([super_target]), weight=1)
            if stoichiometries != None:
                stoichiometries[hid] = {}
                stoichiometries[hid][super_target] = 1.0
                stoichiometries[hid][t] = 1.0

    ## TODO: Convert to Graph
    if stoichiometries == None:
        return super_source,super_target
    else:
        return super_source,super_target,stoichiometries

def median(vals):
    vals.sort()
    if len(vals) % 2 != 0:
        return vals[int(len(vals)/2)]
    else:
        return (vals[len(vals)/2] + vals[len(vals) / 2 + 1])/2

if __name__ == '__main__':
    main(sys.argv)
