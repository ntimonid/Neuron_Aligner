# Initialization of utilized libraries

from cfg import *

# +
# main_path = os.path.abspath('../../')
# lib_dir = os.path.join(main_path,'Libraries')
# data_repository = os.path.join(main_path,'Data Repositories/mouse_connectivity')
# mouselight_dir = os.path.join(main_path,'Data Repositories/Mouselight/json')
# braintell_dir = os.path.join(main_path,'Data Repositories/Braintell')

# sys.path.append(main_path)
# sys.path.append(lib_dir)

# from Libraries.cfg import *
from NeuronMorphology import NeuronMorphology
import rpc_interface
from utils import *


class SBA_interface:

    def __init__(self, sbaHostAlpha = None, sbaHostIncf = None):

        if sbaHostAlpha is None:
            sbaHostAlpha = 'https://neuroinformatics.nl/sba-alpha/www'
        if sbaHostIncf is None:
            sbaHostIncf = 'https://sba-dev.incf.org'

        self.sbaInterface = rpc_interface.RpcInterface(
          'sbaInterface',
          remoteUrl = '{}/{}?{}'.format(sbaHostAlpha,'composer','template=ABA_v3&scene={"background":"FFFFFF"}'),
          interfaceScript = sbaHostIncf+'/js/rpc-interface.js'
        )

        self.soma_thr = None
        self.cpd_params = None

    def set_cpd_params(self,cpd_params = None, soma_thr = None):
        if cpd_params is None:
            self.cpd_params = {'max_it': 3, 'flag_in' : [1,1,-1], 'tol' : 0.001, 'branch_constraint': False}
        else:
            self.cpd_params = cpd_params
        if soma_thr is None:
            self.soma_thr = 30
        else:
            self.soma_thr = soma_thr

    def create_soma_object(self, mySomaLocations, myRegion):
        self.mySomaLocations = mySomaLocations

        markers = []
        for fid,info in mySomaLocations.items():
            coord = info['coord']
            leftHemisphere = coord[2]<11400/20
            markers.append({
                'coord':[1e-2*(coord[0]),1e-2*(coord[1]),1e-2*(coord[2] if leftHemisphere else 1140-coord[2])],
                'color':'#FFFF00' if fid[0] == 'A' else '#FF0000',
                'onclick': {
                    "result": {
                        "db": info['db'],
                        "fname": info['fname'],
                        "flipped": not(leftHemisphere)
                    }
                }
            })

        sbaCommand_center = {
          "method":"Composer.scatter3d",
          "params": {
            "@type": "bas:VectorGraphics",
            "name": u"Soma locations of {} neurons".format(myRegion),
            "space": "sba.ABA_v3(mm,PIR,corner)",
            "style": {
              "marker": {
                "size": 0.05,
              }
            },
            "markers":markers
          }
        }

        return sbaCommand_center



def compare_source_to_targets(source_neuron, target_neuron_ids, cpd_params = None, flip = False):

    if cpd_params is None:
        cpd_params = {'max_it': 30, 'flag_in' : [1,1,-1], 'tol' : 0.001, 'branch_constraint': False}

    # Source pre-processing
    source_morpho_cls = NeuronMorphology(neuronDict = source_neuron)
    source_name = source_neuron['fname'].split('.')[0]
    source_morpho_cls.transform(out_orientation = ['um(10)','PIR','corner'])
    source_minor_lines, source_minor_points = source_morpho_cls.subsampledCopy([1,2], minDistance = 1e10)
    source_soma_point = [line[1] for lineIdx,line in enumerate(source_minor_lines) if line[0] == 1][0]
    X = np.array(source_minor_points)
    X[0,:] = X[1,:]
    X_c = X - X[source_soma_point,:]

    Affinity_dict = OrderedDict(); Affinity_dict_trs = OrderedDict()
    Affinity_dict[source_name] = OrderedDict(); Affinity_dict_trs[source_name] = OrderedDict()

    for target_neuron_id in target_neuron_ids:

        # The following five lines are in case we want to not load stored neurons but download them on the fly instead ...
        db = 'mouselight' if 'AA' in target_neuron_id else 'braintell'
        neuronPath = 'https://neuroinformatics.nl/HBP/neuronsreunited-viewer/{}_json_gz/{}.json.gz'.format(db,target_neuron_id)
        file_content = requests.get(neuronPath).content

        target_neuron = json.loads(zlib.decompress(file_content, 16+zlib.MAX_WBITS))
        target_morpho_cls = NeuronMorphology(neuronDict = target_neuron)
        target_morpho_cls.transform(out_orientation = ['um(10)','PIR','corner'])
        target_name = target_neuron_id #target_morpho_cls.filename.split('.')[0]

        target_minor_lines, target_minor_points = target_morpho_cls.subsampledCopy([1,2], minDistance = 1e10)
        Y = np.array(target_minor_points)
        target_soma_point = [line[1] for lineIdx,line in enumerate(target_minor_lines) if line[0] == 1][0]
        Y[0,:] = Y[1,:]
        if Y[target_soma_point, 2] > 1140//2:  # needs flipping ...
                Y[:,2] = 1140 - Y[:,2]
        Y_c = Y - Y[target_soma_point, :]

        # Time for registration ...
        reg = RigidRegistration(**{'X': X_c, 'Y': Y_c},\
                                tolerance = cpd_params['tol'], max_iterations = cpd_params['max_it'],
                                flag_in = cpd_params['flag_in'])
        TY, (s_reg, R_reg, t_reg) = reg.register()

        match_1_to_2 = np.argmax(reg.P,axis = 1)
        X_c_remap = X_c[match_1_to_2,:]
        MSE_trs = np.linalg.norm(X_c_remap[:,0:3] - TY[:,0:3])
        X_remap = X[match_1_to_2,:]
        MSE = np.linalg.norm(X_remap[:,0:3] - Y[:,0:3])

        Affinity_dict[source_name][target_name] = MSE_trs # MSE

    return Affinity_dict


### Call the click handler
def neuron_finder(result, sba_interface, mySomaLocations, neuriteLengthDistribution):

    sbaInterface = sba_interface.sbaInterface
    nld_list = neuriteLengthDistribution.keys()

    mySomaLocations_cpy = deepcopy(mySomaLocations)
    for key,val in mySomaLocations_cpy.items():
        soma_cord = np.array(mySomaLocations_cpy[key]['coord'])
        if soma_cord[2] > 1140//2:
            mySomaLocations_cpy[key]['coord'][2] = 1140 - soma_cord[2]

    db = result['db']
    neuronName = result['fname']
    infile = neuronName + '.json.gz'
    if db == 'mouselight':
        orient = 'LIP'; scale = 'um'; position = 'corner'; encoding = "gzip"
#         neuronFile = os.path.join(mouselight_dir,infile)
    elif db == 'braintell':
        orient = 'PIR'; scale = 'um'; position = 'corner'; encoding = "gzip"
#         neuronFile = os.path.join(braintell_dir,infile)

    neuronPath = 'https://neuroinformatics.nl/HBP/neuronsreunited-viewer/{}_json_gz/{}.json.gz'.format(db,neuronName)
    file_content = requests.get(neuronPath).content
    flipped = ''
    if result['flipped']:
        file_content,encoding, flipped = flip(file_content, orient)
        source_neuron_dict = json.loads(file_content)
    else:
        source_neuron_dict = json.loads(zlib.decompress(file_content, 16+zlib.MAX_WBITS))
    source_neuron_dict['fname'] = neuronName
    color = '#00FF00'

    # Here I shall try to mimick compare neurons ...
    sbaCommand_click = {
      "method":"Composer.message",
      "params" : {
        "message":"You clicked on an object with params {} {} {} {} {}".format(neuronName, neuronPath, orient,
                                                                               scale, position)
      }
    }

    sbaCommand_load_morpho = {
          "method":"Composer.import",
          "params": {
            "@type": "bas:morphology",
            "name": "{}{}.json".format(neuronName,flipped),
            "contents": file_content,
            "mediatype": "application/vnd.hbp.movi+json",
            "encoding": encoding,
            "space": "sba:ABA_v3({},{},{})".format(scale,orient,position),
            "style": {
              "axonColor": '{}'.format(color)
            }
          }
    }

    sbaInterface.send(sbaCommand_click)
    sbaInterface.send(sbaCommand_load_morpho)

    if sba_interface.soma_thr is None and sba_interface.cpd_params is None:  # need to intialize ...
        sba_interface.set_cpd_params()
    soma_thr, cpd_params = sba_interface.soma_thr, sba_interface.cpd_params

    neuronName_abbrev = braintell_2_nld(nld_list, neuronName)
    source_cords = mySomaLocations_cpy[neuronName_abbrev]['coord']
    search_list = []

    for neuron_info in mySomaLocations_cpy.values():
        target_cords,neuron_id = np.array(neuron_info['coord']), neuron_info['fname']
        if neuron_id == neuronName: continue
        soma_distance = np.linalg.norm(source_cords - target_cords)
        if soma_distance >= sba_interface.soma_thr: continue
        search_list.append(neuron_id)

    sbaCommand_wait = {
      "method":"Composer.message",
      "params" : {
        "message":"Searching for the most similar neuron ..."
      }
    }
    sbaInterface.send(sbaCommand_wait)

    Affinity_dict = compare_source_to_targets(source_neuron_dict, search_list, cpd_params)

    neuron_1 = list(Affinity_dict.keys())[0]
    trg_neurons = list(Affinity_dict[neuron_1].keys())
    mses = [val for key,val in Affinity_dict[neuron_1].items()] #['MSE']

    if len(mses) == 0:
        nearest_neuron = "None"
        sbaCommand_nearest_click = {
          "method":"Composer.message",
          "params" : {
            "message":"No proximals neurons were found. Increase search radius"
          }
        }
        sbaInterface.send(sbaCommand_nearest_click)
    else:
        nearest_neuron = trg_neurons[numpy.argmin(mses)]

        sbaCommand_nearest_click = {
          "method":"Composer.message",
          "params" : {
            "message":"The nearest neuron is: {}".format(nearest_neuron)
          }
        }
        sbaInterface.send(sbaCommand_nearest_click)

        file_content, orient, scale, position, encoding, hemisphere = GetNeuronInfo(nearest_neuron)

        flipped = ''

        if hemisphere == 'right':
            file_content,encoding, flipped = flip(file_content, orient)

        sbaCommand_nearest_morpho = {
                  "method":"Composer.import",
                  "params": {
                    "@type": "bas:morphology",
                    "name": "{}{}.json".format(nearest_neuron,flipped),
                    "contents": file_content,
                    "mediatype": "application/vnd.hbp.movi+json",
                    "encoding": encoding,
                    "space": "sba:ABA_v3({},{},{})".format(scale,orient,position),
                    "style": {
                      "axonColor": '{}'.format(color)
                    }
                  }
            }
        sbaInterface.send(sbaCommand_nearest_morpho)



# def make_population(source_neuron, target_neuron_ids, cpd_params = None, flip = False):
#
#
#     if cpd_params is None:
#         cpd_params = {'max_it': 30, 'flag_in' : [1,1,-1], 'tol' : 0.001, 'branch_constraint': False}
#
#     # Source pre-processing
#     source_morpho_cls = NeuronMorphology(neuronDict = source_neuron)
#     source_name = source_neuron['fname'].split('.')[0]
#     source_morpho_cls.transform(out_orientation = ['um(10)','PIR','corner'])
#     source_minor_lines, source_minor_points = source_morpho_cls.subsampledCopy([1,2], minDistance = 1e10)
#     source_soma_point = [line[1] for lineIdx,line in enumerate(source_minor_lines) if line[0] == 1][0]
#     X = np.array(source_minor_points)
#     X[0,:] = X[1,:]
#     X_c = X - X[source_soma_point,:]
#
#     Affinity_dict = OrderedDict(); Affinity_dict_trs = OrderedDict()
#     Affinity_dict[source_name] = OrderedDict(); Affinity_dict_trs[source_name] = OrderedDict()
#
#     for target_neuron_id in target_neuron_ids:
#
#         # Target pre-processing
#         target_morpho_cls = NeuronMorphology(neuronFile = target_neuron_id)
#         target_morpho_cls.transform(out_orientation = ['um(10)','PIR','corner'])
#         target_name = target_morpho_cls.filename.split('.')[0]
#         target_minor_lines, target_minor_points = target_morpho_cls.subsampledCopy([1,2], minDistance = 1e10)
#         Y = np.array(target_minor_points)
#         target_soma_point = [line[1] for lineIdx,line in enumerate(target_minor_lines) if line[0] == 1][0]
#         Y[0,:] = Y[1,:]
#         if Y[target_soma_point, 2] > 1140//2:  # needs flipping ...
#                 Y[:,2] = 1140 - Y[:,2]
#         Y_c = Y - Y[target_soma_point, :]
#
#         # Time for registration ...
#         reg = RigidRegistration(**{'X': X_c, 'Y': Y_c},\
#                                 tolerance = cpd_params['tol'], max_iterations = cpd_params['max_it'],
#                                 flag_in = cpd_params['flag_in'])
#         TY, (s_reg, R_reg, t_reg) = reg.register()
#
#         match_1_to_2 = np.argmax(reg.P,axis = 1)
#         X_c_remap = X_c[match_1_to_2,:]
#         MSE_trs = np.linalg.norm(X_c_remap[:,0:3] - TY[:,0:3])
#         X_remap = X[match_1_to_2,:]
#         MSE = np.linalg.norm(X_remap[:,0:3] - Y[:,0:3])
#
#         Affinity_dict[source_name][target_name] = MSE_trs # MSE
#
#     return Affinity_dict
