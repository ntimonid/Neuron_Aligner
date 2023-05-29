from cfg import *

def load_useful_variables(data_repository):

    try:
        with open(os.path.join(data_repository, 'ccf3_acr2id.json')) as fp:
            acr2id = json.load(fp)
        with open(os.path.join(data_repository,'ancestorsById.json')) as fp:
            ancestorsById = json.load(fp)
    except:
        from allensdk.core.reference_space_cache import ReferenceSpaceCache

        reference_space_key = 'annotation/ccf_2017'
        rspc = ReferenceSpaceCache(10, reference_space_key, manifest='manifest.json')
        tree = rspc.get_structure_tree(structure_graph_id=1)
        acr2id = tree.get_id_acronym_map()
        ancestorsById = tree.get_ancestor_id_map()

    id2acr = { id:acr for acr,id in acr2id.items() }
    if 0 not in id2acr:
      acr2id['[background]'] = 0
      id2acr[0] = '[background]'

    with open(os.path.join(data_repository,'acr_to_morpho_id_new.pkl'), 'rb') as infile:
        acr_to_morpho_id = pk.load(infile)

    neuriteLengthDistribution = {}
    databases = ['braintell','mouselight']
    for dbName in databases:
          with open(os.path.join(data_repository, 'neuriteLengthDistribution({}).json'.format(dbName))) as fp:
                neuriteLengthDistribution.update(json.load(fp))

    return acr2id, ancestorsById, neuriteLengthDistribution, acr_to_morpho_id



def braintell_2_nld(nld_list, morpho_id):
    if 'AA' in morpho_id:
        return morpho_id
    else:
        try:
            tmp = morpho_id.split('_reg')[0]
            matches =  [key for key in nld_list if tmp.split('_')[2] in key]
            return matches[0]
        except:
            return -1


def _micron_checker(somaCoord, orient = 'PIR'):
      if somaCoord[0] < 10 and somaCoord[1] < 10 and somaCoord[2] < 10:
          res = 'mm'
          unit = 1
      else:
          res = 'um'
          if orient == 'PIR':
              P,I,R = (528, 320, 456)
              unit = 1 if (somaCoord[0]>P or somaCoord[1]>I or somaCoord[2]>R) else 25
          elif orient == 'LIP':
              L,I,P = (450, 320, 528)
              unit = 1 if (somaCoord[0]>L or somaCoord[1]>I or somaCoord[2]>P) else 25
      return res,unit

def GetNeuronInfo(neuronName):

        braintell_link = 'https://neuroinformatics.nl/HBP/braintell-viewer/fname2soma.bas(sba.ABA_v3(mm,RAS,ac)).json'
        response = urllib.request.urlopen(braintell_link)
        file_content = response.read()
        out = json.loads(file_content)

        if 'AA' in neuronName:
            mode = 'mouselight'
            orient = 'LIP'; position = 'corner'; encoding = "gzip"
            # scale = 'um';
            infile = neuronName +'.json.gz'
        else:
            mode = 'braintell'
            orient = 'PIR'; position = 'corner'; encoding = "gzip"
            # scale = 'um';
            neuronName = [key for key in out.keys() if neuronName in key][0]
            infile = neuronName + '.json.gz'

        neuronPath = 'https://neuroinformatics.nl/HBP/neuronsreunited-viewer/{}_json_gz/{}.json.gz'.format(mode,neuronName)
        contents = requests.get(neuronPath).content #b64encode(requests.get(neuronPath).content).decode('utf-8')
        neuron_dict = json.loads(zlib.decompress(contents, 16+zlib.MAX_WBITS))
        somaIdx = [line[1] for lineIdx,line in enumerate(neuron_dict['treeLines']['data']) if line[0] == 1][0]
        somaCoord = neuron_dict['treePoints']['data'][somaIdx]

        res,unit = _micron_checker(somaCoord, orient)
        if unit == 25:
            scale = '{}({})'.format(res,unit)
        else:
            scale = '{}'.format(res)

        far_right = [11400 if unit == 1 else 11400//25][0]
         
        if orient == 'PIR':
            if somaCoord[2] > far_right//2: # time to flip
                hemisphere = 'right'
            else:
                hemisphere = 'left'
        elif orient == 'LIP':
            if somaCoord[0] < far_right//2: # time to flip
                hemisphere = 'right'
            else:
                hemisphere = 'left'

        return contents, orient, scale, position, encoding, hemisphere


def get_soma_locations(myRegion, neuriteLengthDistribution, acr2id, ancestorsByid):

    mySomaLocations = {}

    un_num = lambda x : x.split('1')[0].split('2/3')[0].split('4')[0].split('5')[0].split('6a')[0].split('6b')[0]
    cortexId = acr2id['Isocortex']
    myId = acr2id[myRegion]

    for fid,nld in neuriteLengthDistribution.items():

        if 'AA' in fid:
            db = 'mouselight'
        else:
            db = 'braintell'
        somaRegion = nld['soma']['region']
        if somaRegion == '[?]' or un_num(nld['soma']['region']) != un_num(nld['soma']['correctedRegion']): continue
        id = acr2id[somaRegion] if somaRegion in acr2id else 0
        ancestors = ancestorsByid[str(id)] if str(id) in ancestorsByid else []
        if myId in ancestors: #== id: #
            mySomaLocations[fid] = {
                'coord': nld['soma']["coord(um(10),PIR,corner)"],
                'fname': nld['fname'],
                'db': db
            }

    return mySomaLocations


def flip(file_content, orient):

    flipped = '(flipped)'
    neuron = json.loads(zlib.decompress(file_content, 16+zlib.MAX_WBITS) )
    lrPos = 1 + orient.find('L') + orient.find('R')

    for coord in neuron['treePoints']['data']:
        coord[lrPos] = 11400-coord[lrPos]

    file_content = json.dumps(neuron)
    encoding = ""

    return file_content, encoding, flipped




def CallSBA(neuronName, neuronPath, in_orientation = None, mode = 'mouselight',
            color = '#00FF00', flipped = False, interface = True):

        # main function for calling SBA using python inputs
        # neuronName: the name of the imported morphology or streamline, as you want it to show up at the SBA website
        # neuronPath: here you can either a) give the url of a morphology from a repository of interest, or
        #             directly give the morphology you are currently analysing as input in a python-dictionary-style format
        if type(neuronPath) == dict:
            if in_orientation is None:
                print('please specify the orientation system of your input data')
                return -1
            out_orientation = ['mm','RAS','ac']
            encoding = ""
            scale,orient,position = out_orientation
            Q = CAS.convertAllenSpace(in_orientation,out_orientation)
            points = np.array(neuronPath['treePoints']['data'])
            neuronPath['treePoints']['data'] = (np.matmul(Q[0:3,0:3],points[:,0:3].T).T +  Q[0:3,3]).tolist()
            file_content = json.dumps(neuronPath)

        elif mode.lower() == 'mouselight':
            orient = 'LIP'; scale = 'um'; position = 'corner'; encoding = "gzip"
            file_content = requests.get(neuronPath).content
        elif mode.lower() == 'braintell':
            orient = 'PIR'; scale = 'um'; position = 'corner'; encoding = "gzip"
            file_content = requests.get(neuronPath).content
        elif mode.lower() == 'streamline':
            streamline_x3d = streamlines2x3d(neuronPath)
        elif mode.lower() == 'local':
            orient = 'RAS'; scale = 'mm'; position = 'ac'; encoding = ""
            with open(neuronPath,'rt') as fp: file_content = fp.read()

        if flipped is False:
            flipped = ''
        else:
            file_content, encoding, flipped = flip(file_content, orient)

        if mode.lower() == 'streamline':
            sbaCommand_load_morpho = {
                "method": "Composer.import",
                "params": {
                    "message": "",
                    "mime": 'model/x3d',
                    "name": 'streamlines_' + str(neuronName) + '.bas{sba.ABA_v3^corner[PIR,um]}.x3d',
                    "contents": streamline_x3d,
                }
            }
        else:
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

        if interface is True:
            # latest addition to include two versions of the CallSBA: One where you just create the command,
            # and one where you create a new interface and load the command
#             sba_inter = SBA_interface()
            sba_inter.sbaInterface.send(sbaCommand_load_morpho)

        return sbaCommand_load_morpho

# ## DownloadAxons function
# This is the inbuilt function used inside the tree2volume and JSON_the_cluster functions in order to download either the mouselight neuronal morphology reconstructions or the Allen Mouse Brain Connectivity Atlas-based streamlines.
# ### input:
# experimentId = Id of the experiment corresponding to the morphology/streamline of interest.
# mode = 'mouselight', 'braintell, ''streamlines', or 'local (Clasca neurons from the Consortium)' for the respective data modality. Default is 'mouselight'
# rgb = the colormap to be used for visualizing the imported morphologies or streamlines. Default is '#FFBB00'

def DownloadAxons(experimentId = None, mode = 'mouselight', rgb = '#FFBB00'): # Needs an update

    if mode == 'streamlines':
        save_path = os.path.join(main_path, 'Backup_Code/streamlines_tmp')
        infile = 'streamlines_{}.json.gz'.format(experimentId)
        streamline_link = 'https://neuroinformatics.nl/HBP/allen-connectivity-viewer/json/{}'.format(infile)
    elif mode == 'mouselight':
        save_path = os.path.join(main_path, 'Data Repositories/Mouselight/json')
        if type(experimentId) == int: # id not given directly
            experimentId += 1
            infile = 'AA{:04d}'.format(experimentId)
            experimentId = infile
        if '.json' not in experimentId:
            infile = '{}.json.gz'.format(experimentId)
        else:
            infile = experimentId
        streamline_link = 'https://neuroinformatics.nl/HBP/neuronsreunited-viewer/{}_json_gz/{}'.format(mode,infile)
    elif mode == 'braintell':
        save_path = os.path.join(main_path, 'Data Repositories/Braintell')
        braintell_link = 'https://neuroinformatics.nl/HBP/braintell-viewer/fname2soma.bas(sba.ABA_v3(mm,RAS,ac)).json'
        response = urllib.request.urlopen(braintell_link)
        file_content = response.read()
        out = json.loads(file_content)
        if type(experimentId) == int: # id not given directly
            cnt = 0
            for key in out.keys():
                if cnt == experimentId:
                    experimentId = key
                    break
                cnt+=1
            infile = '{}.json.gz'.format(experimentId)
        if '.json.gz' not in experimentId:
            infile = '{}.json.gz'.format(experimentId)
        else:
            infile = experimentId
        streamline_link = 'https://neuroinformatics.nl/HBP/neuronsreunited-viewer/{}_json_gz/{}'.format(mode,infile)

    infile2 = os.path.join(save_path,infile)
    print(infile, infile2, streamline_link)
    try:
        if os.path.exists(infile2) is False:
            wget.download(streamline_link,infile2)
        if '.json.gz' in infile2:
            with gzip.open(infile2, 'rb') as fp:
                file_content = fp.read()
        elif '.json' in infile2:
            with open(infile2, 'r') as fp:
                file_content = fp.read()
        out = json.loads(file_content)
    except:
        out = -1

    return out


def Acro2Morpho(acronym, mode = 'mouselight', morpho_dict = None):

    # For this function you will need to download the allensdk library: https://allensdk.readthedocs.io/en/latest/
    from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache

    # This function will allow us to link brain area acronyms, as denoted by the ARA notation,  to specific experiment Ids

    mcc = MouseConnectivityCache()
    structure_tree = mcc.get_structure_tree()

    if mode == 'streamlines':
        acro_id = structure_tree.get_structures_by_acronym([acronym])[0]
        experiment_list = mcc.get_experiments(injection_structure_ids=[acro_id['id']])
        experiment_ids = [val['id'] for val in experiment_list]
    else:
        if morpho_dict is None:
            print('Error! Give me a dict of morphology files')
            return -1
        experiment_ids = [val for val in morpho_dict[acronym]]

    return experiment_ids


def is_positive_semi_definite(R):
    if not isinstance(R, (np.ndarray, np.generic)):
        raise ValueError('Encountered an error while checking if the matrix is positive semi definite. \
            Expected a numpy array, instead got : {}'.format(R))
    return np.all(np.linalg.eigvals(R) > 0)
