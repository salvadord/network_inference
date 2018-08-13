# infer netpyne model conns usign graph-tool
import json
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

# infer high-level connectivity using graph-tool
def infer(net, outputFile):

    # iterate over conns and create graph
    from graph_tool import all as gt

    cells,conns = [], []
    g = gt.Graph()
    pops = g.new_vertex_property("string")
    pos = g.new_vertex_property("vector<double>")

    for cell in net['cells']:
        cells.append(g.add_vertex())
        pops[cells[-1]] = cell['tags']['pop']
        pos[cells[-1]] = [cell['tags'][p] for p in ['x', 'y', 'z']]

    for cell in net['cells']:
        postGid = cell['gid']
        for conn in cell['conns']:
            preGid = conn['preGid']
            if preGid != 'NetStim':
                conns.append(g.add_edge(cells[int(preGid)], cells[int(postGid)]))


    gt.graph_draw(g, vertex_text=g.vertex_index, pos=pos, vertex_font_size=5,output_size=(1000, 1000), output=outputFile+'_graph.png')

    #bv, be = betweenness(g)

    state = gt.minimize_blockmodel_dl(g)
    state.draw(pos=pos, output=outputFile+'_state.png')
    e = state.get_matrix()
    print(e)
    plt.matshow(e.todense().T)
    plt.savefig(outputFile+'_mat.png')


if __name__ == '__main__':
    # params
    model = 'net1'
    inputFile = model+'.json'
    outputFile = model+'_infer'

    #runModel(model)
    net=loadData(inputFile)
    infer(net, outputFile)



