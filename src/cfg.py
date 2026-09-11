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
import copy
from copy import deepcopy
from collections import OrderedDict

try:
    import pandas as pd
except ImportError:
    pd = None
import numpy as np
import pickle as pk

try:
    import ipywidgets as widgets
    from IPython.core.display import display, HTML
    from IPython.display import Javascript, clear_output
except ImportError:
    class Dummy:
        def __getattr__(self, *args, **kwargs): return self
        def __call__(self, *args, **kwargs): return None
    widgets = Dummy()
    display = lambda *a, **kw: None
    HTML = lambda *a, **kw: None
    Javascript = lambda *a, **kw: None
    clear_output = lambda *a, **kw: None
from json import dumps as json_encode, loads as json_decode
from base64 import b64encode,b64decode
try:
    from urlparse import urljoin  # Python2
except ImportError:
    from urllib.parse import urljoin

