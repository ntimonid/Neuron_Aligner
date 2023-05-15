from __future__ import unicode_literals

import os
import sys
import urllib
import wget
import traceback
import uuid
import gzip
import zlib
import json
import numpy
import requests
import re
import nrrd
import copy

import ipywidgets as widgets
import pandas as pd
import numpy as np
import pickle as pk

from base64 import b64encode
from json import dumps as json_encode
from IPython.core.display import display, HTML
from copy import deepcopy
from collections import OrderedDict

from IPython.display import Javascript, clear_output #display,
from json import dumps as json_encode, loads as json_decode
from base64 import b64encode,b64decode
try:
    from urlparse import urljoin  # Python2
except ImportError:
    from urllib.parse import urljoin

from cpd_registration import RigidRegistration
