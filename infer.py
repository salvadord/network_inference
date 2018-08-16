# infer netpyne model conns usign graph-tool
import matplotlib; matplotlib.use('Agg')

import json
from graph_tool import all as gt
from time import time


import matplotlib.pyplot as plt

# run model
def runModel(model):
    import importlib
    importlib.import_module(model)

# load model output data()
def loadData(filename):
    with open(filename, 'r') as f: 
        data = json.load(f)

    net = data['net']
    return net


if __name__ == '__main__':
    start = time()

    # --------------------------------------------------------------------------
    # params
    # --------------------------------------------------------------------------
    model = 'net2'
    inputFile = 'data/'+model+'.json'
    outputFile = 'data/'+model+'_infer'

    # --------------------------------------------------------------------------
    # run model and load data
    # --------------------------------------------------------------------------
    #runModel(model)
    net=loadData(inputFile)  

    # --------------------------------------------------------------------------
    # infer connectivity
    # --------------------------------------------------------------------------
    from netpyne import analysis
    colors = analysis.colorList
    colorDict = {pop: i for i,pop in enumerate(net['pops'])} #(color+[1.0] for pop,color in zip(net['pops'],colors)}

    # iterate over conns and create graph
    cells, conns = [], []
    g = gt.Graph()
    pops = g.new_vertex_property("string")
    pop_ids = g.new_vertex_property("int")
    locs = g.new_vertex_property("vector<double>")
    weights = g.new_edge_property("double")
    g.vertex_properties['pops'] = pops
    g.vertex_properties['pops_ids'] = pop_ids
    g.vertex_properties['locs'] = locs
    g.edge_properties['weights'] = weights

    for origCell in net['cells']:
        if origCell['tags']['pop'] not in ['bkg']:
            cells.append(g.add_vertex())
            pop = origCell['tags']['pop']
            pop_id = colorDict[pop]
            pops[cells[-1]] = pop
            pop_ids[cells[-1]] = pop_id
            locs[cells[-1]] = [origCell['tags'][p] for p in ['x', 'y', 'z']]

    for origCell in net['cells']:
        postGid = origCell['gid']
        for origConn in origCell['conns']:
            preGid = origConn['preGid']
            if preGid != 'NetStim':
                conns.append(g.add_edge(cells[int(preGid)], cells[int(postGid)]))
                weights[conns[-1]] = origConn['weight']

    # plot graph (like plot2Dnet) # pos=g.vp.locs 
    #pos = gt.planar_layout(g)  #radial_tree_layout(g, 0) #arf_layout(g) # fruchterman_reingold_layout(g) 
    gt.graph_draw(g, pos=g.vp.locs,
                    vertex_text=g.vertex_index, 
                    vertex_size=5, 
                    vertex_fill_color=pop_ids, 
                    vcmap=matplotlib.cm.jet, 
                    vcnorm=1,
                    edge_color=g.ep.weights, 
                    ecmap=(matplotlib.cm.jet, .6),
                    edge_pen_width=gt.prop_to_size(g.ep.weights, 0.05, 0.1, power=1), 
                    vertex_font_size=5, 
                    output_size=(100, 1000), 
                    output=outputFile+'_graph.pdf')

    # caculate betweeness
    # bv, be = betweenness(g)

    # plot states (cell groupings)
    #state = gt.minimize_blockmodel_dl(g)
    #state = gt.minimize_nested_blockmodel_dl(g)
    #state_ndc = gt.minimize_nested_blockmodel_dl(g, deg_corr=False)
    #state_dc  = gt.minimize_nested_blockmodel_dl(g, deg_corr=True)
    #state_dc  = gt.minimize_nested_blockmodel_dl(g, deg_corr=True, mcmc_equilibrate_args={'force_niter':100, 'mcmc_args':dict(niter=5)})
    state_dc_w  = gt.minimize_nested_blockmodel_dl(g, deg_corr=True)#, state_args=dict(recs=[g.ep.weights], rec_types=['real-exponential']))

    # print("Non-degree-corrected DL:\t", state_ndc.entropy())
    # print("Degree-corrected DL:\t", state_dc.entropy())
    # print(u"ln \u039b: ", state_dc.entropy() - state_ndc.entropy())

    state=state_dc_w

    state.print_summary()
    state.draw(pos=g.vp.locs, output=outputFile+'_state.png')  # vertex_shape="pie", vertex_pie_fractions=pv, vertex_pie_colors --> use to represent pops? no
    
    # analysis
    levels = state.get_levels()
    groupCells = {}
    groupPops = {}
    for ilevel, level in enumerate(levels):
        print(level)
        b = level.get_blocks()
        numGroups = len(list(state.levels[0].wr))
        groupCells[ilevel] = {}
        groupPops[ilevel] = {}

        for group in range(numGroups):
            cells = [i for i, x in enumerate(b) if x == group]
            groupCells[ilevel][group] = cells
            if len(cells) > 0:
                groupPops[ilevel][group] = {origPopName: len(list(set(cells) & set(origPop['cellGids'])))/len(cells) for origPopName,origPop in net['pops'].items()}

        # bg = level.get_bg()
        # ers = level.mrs    # edge counts
        # print(ers.a)
        # nr = level.wr      # node counts
        # print(nr.a)


        # plot conn matrix
        e = level.get_matrix()
        connMat = e.todense().T
        plt.matshow(connMat)  
        plt.savefig(outputFile+'_level_'+str(ilevel)+'_mat.png')

        # TODO: ADD LABELS TO CONN MATRIX AXES (G1, G2, ...)

        # TODO: PIE CHARTS FOR EACH GROUP SHOWING % OF EACH POP

        # GRAPH WITH LOCATIONS OF CELLS AND PIE CHART SHOWING GROUP COLOR + POP COLOR (SAME COLOR IF SAME)

    # plot elapsed time
    end = time()
    print('elapsed time:%d s'%(end - start))

