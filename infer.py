# infer netpyne model conns usign graph-tool
import matplotlib; matplotlib.use('Agg')

import numpy as np
import json
from graph_tool import all as gt
from time import time
from PIL import Image, ImageFont, ImageDraw

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
    model = 'm1'
    inputFile = 'data/'+model+'.json'
    outputFile = 'data/'+model+'_infer'
    loadState = True

    # --------------------------------------------------------------------------
    # run model and load data
    # --------------------------------------------------------------------------
    print('Running and/or loading network model...')
    #runModel(model)
    net=loadData(inputFile)  

    # --------------------------------------------------------------------------
    # convert network to graph
    # --------------------------------------------------------------------------
    from netpyne import analysis
    colors = analysis.colorList
    colorDict = {pop: i for i,pop in enumerate(net['pops'])} #(color+[1.0] for pop,color in zip(net['pops'],colors)}

    # iterate over conns and create graph
    print('  elapsed time: %d s'%(time() - start))
    print('Converting loaded network to graph...')
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

    preGidKey = 0 #'preGid' # 0
    weightKey = -2 #'weight' # -2
    for origCell in net['cells']:
        postGid = origCell['gid']
        for origConn in origCell['conns']:
            preGid = origConn[preGidKey]
            if preGid != 'NetStim':
                conns.append(g.add_edge(cells[int(preGid)], cells[int(postGid)]))
                weights[conns[-1]] = origConn[weightKey]

    # plot graph (like plot2Dnet) # pos=g.vp.locs 
    #pos = gt.planar_layout(g)  #radial_tree_layout(g, 0) #arf_layout(g) # fruchterman_reingold_layout(g) 
    plotGraph = 0
    if plotGraph:
        print('Plotting node and edges graph...')
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

    print('  elapsed time: %d s'%(time() - start))


    # --------------------------------------------------------------------------
    # run inference SBM algorithm to extract groupings
    # --------------------------------------------------------------------------
    if loadState:
        print('Loading the nested stochastic block model...')
        import pickle
        with open('data/m1_500k_state.pkl', 'rb') as f:
            state = pickle.load(f)['state']
    else:
        print('Fitting the nested stochastic block model, by minimizing its description length...')
        #state = gt.minimize_blockmodel_dl(g)
        #state = gt.minimize_nested_blockmodel_dl(g)
        #state_ndc = gt.minimize_nested_blockmodel_dl(g, deg_corr=False)
        #state_dc  = gt.minimize_nested_blockmodel_dl(g, deg_corr=True)
        #state_dc  = gt.minimize_nested_blockmodel_dl(g, deg_corr=True, mcmc_equilibrate_args={'force_niter':100, 'mcmc_args':dict(niter=5)})
        state_dc_w  = gt.minimize_nested_blockmodel_dl(g, deg_corr=True, state_args=dict(recs=[g.ep.weights], rec_types=['real-exponential']))

        # print("Non-degree-corrected DL:\t", state_ndc.entropy())
        # print("Degree-corrected DL:\t", state_dc.entropy())
        # print(u"ln \u039b: ", state_dc.entropy() - state_ndc.entropy())


        # state analysis
        print('  elapsed time: %d s'%(time() - start))
        print('Analysing states of nested block model ...')
        state=state_dc_w

    state.print_summary()

    # plot nested SBM
    pos,t,tpos = state.draw(pos=g.vp.locs, output=outputFile+'_state.png', output_size=(1200, 1200)) #, layout ="sfdp")  # vertex_shape="pie", vertex_pie_fractions=pv, vertex_pie_colors --> use to represent pops? no
 

    # --------------------------------------------------------------------------
    # analyze groupings and conn matrix in each level
    # --------------------------------------------------------------------------
    levels = state.get_levels()
    groupCells = {}
    groupPops = {}
    groupOrder = {}
    for ilevel, level in enumerate(levels):
        print(level)
        b = level.get_blocks()
        numGroups = len(level.wr.a)
        groupCells[ilevel] = {}
        groupPops[ilevel] = {}
        groupOrder[ilevel] = []

        for group in range(numGroups):
            cells = [i for i, x in enumerate(b) if x == group]
            groupCells[ilevel][group] = cells
            if len(cells) > 0:
                groupPops[ilevel][group] = {origPopName: len(list(set(cells) & set(origPop['cellGids'])))/len(cells) for origPopName,origPop in net['pops'].items()}

            v=list(groupPops[ilevel][group].values())
            groupOrder[ilevel].append(v.index(max(v)))
        print(groupOrder)

        # plot conn matrix
        e = level.get_matrix()
        if ilevel == 0 and numGroups == len(net['pops']):
            connMat = e.todense()[groupOrder[ilevel]]
            connMat = connMat.T[groupOrder[ilevel]]
            plt.matshow(connMat[groupOrder[ilevel]][groupOrder[ilevel]])  
        else:
            connMat = e.todense().T
            plt.matshow(connMat)

        labels = ['G'+str(i) for i in range(numGroups)]
        plt.xticks(range(numGroups), labels)
        plt.yticks(range(numGroups), labels)

        plt.savefig(outputFile+'_level_'+str(ilevel)+'_mat.png')

        # pie charts for each group
        if ilevel == 0:
            piesPerRow = 6
            sizePerPie = 6
            fig=plt.figure(figsize=(sizePerPie*(piesPerRow), sizePerPie*((numGroups/piesPerRow)+1)))    
            fig.subplots_adjust(right=0.98, left=0.02, top=0.95, bottom=0.05) # Less space on right
            ax = plt.axes([0.05, 0.05, 0.95, 0.95])       
            for igroup, groupPop in groupPops[ilevel].items():
                if numGroups == len(net['pops']):
                    igroupNew = groupOrder[ilevel][igroup]
                else: 
                    igroupNew = igroup
                plt.subplot(int(numGroups/piesPerRow)+1, piesPerRow, igroupNew+1)

                # plot pie chart
                popLabels = [pl for pl in groupPop.keys() if groupPop[pl]>0]
                popFracs = groupPop.values()
                from pylab import mpl
                mpl.rcParams['font.size'] = 20.0
                mpl.rcParams['font.weight'] = 'bold'
                popColors = [colors[colorDict[pl]] for pl in popLabels]
                fracs = [pf for pf in popFracs if pf>0]           
                plt.pie(fracs, labels=popLabels, autopct='%1.0f%%', pctdistance=0.8, labeldistance=1.1, shadow=True, startangle=0, colors=popColors)
                plt.title('G'+str(igroupNew)+' (%d cells)'%(level.wr[igroup]))

            plt.savefig(outputFile+'_level_'+str(ilevel)+'_groupPiecharts.png')


    # --------------------------------------------------------------------------
    # Add labels to nested SBM figure
    # --------------------------------------------------------------------------
    img = Image.open("data/m1_infer_state.png")
    draw = ImageDraw.Draw(img)
    # font = ImageFont.truetype("sans-serif.ttf", 16)
    font = ImageFont.truetype('/Library/Fonts/Arial.ttf', 22)

    numCells = state.levels[0].wr.a.sum()
    numGroups = state.levels[1].wr.a.sum()
    dpi = 100
    baseOffset = 6.0
    factorX=0.1
    factorY=0.1
    popXOffset = {16:-3, 13:-6, 11:-7.5, 24:-4.5, 25:-3.5, 23:-1.5, 10:0, 22:0, 21:0} 
    popYOffset = {11:5.5, 24:3, 25:1, 23:0, 10:4, 22:0, 21:-1} 
    for i in range(numGroups):
        inode = numCells+i
        nodePos = tpos[inode]
        #print(nodePos)
        #print(((nodePos[0]+offset)*dpi,  (nodePos[1]+offset)*dpi))
        gpops = [k for k,v in groupPops[0][i].items() if v > 0] 
        gpops = ','.join(gpops)
        offsetX=baseOffset+popXOffset[i]*factorX if i in popXOffset else baseOffset
        offsetY=baseOffset+popYOffset[i]*factorY if i in popYOffset else baseOffset
        draw.text(((nodePos[0]+offsetX)*dpi, (nodePos[1]+offsetY)*dpi), 'G%d (%s)'%(i,gpops), (0,0,0), font=font)
        #draw.text(((nodePos[0]+offset)*dpi, (nodePos[1]+offset)*dpi), '\n%s'%(gpops), (0,0,0), font=font)
        
    img.save('data/m1_infer_state_labeled.png') 


    #from IPython import embed; embed()

    # GRAPH WITH LOCATIONS OF CELLS AND PIE CHART SHOWING GROUP COLOR + POP COLOR (SAME COLOR IF SAMEm)

    # plot elapsed time
    print('elapsed time: %d s'%(time() - start))

