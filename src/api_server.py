import http.server
import json
import socketserver
import logging
from src import config
from src import sba_api
from legacy.rpc_interface import RpcInterface

logger = logging.getLogger(__name__)

# Simple configuration
PORT = config.PORT

# State
sba_interface = None
neuriteLengthDistribution = None
mySomaLocations = None
dummy_rpc = RpcInterface('dummy', 'url', 'script')

def set_api_state(sba_int, nld, soma_locs):
    global sba_interface, neuriteLengthDistribution, mySomaLocations
    sba_interface = sba_int
    neuriteLengthDistribution = nld
    mySomaLocations = soma_locs

class APIHandler(http.server.BaseHTTPRequestHandler):
    def do_POST(self):
        content_length = int(self.headers['Content-Length'])
        post_data = self.rfile.read(content_length)
        data = json.loads(post_data.decode('utf-8'))
        
        # Handle callback
        if self.path == '/callback':
            packedResponse = data.get('packedResponse')
            callbackName = data.get('callback')
            
            # Unpack response
            response = dummy_rpc.unpackResponse(packedResponse)
            
            # Dispatch to handler
            if callbackName == 'sbaOnclickHandler':
                sba_api.neuron_finder(response.get('result'), sba_interface, mySomaLocations, neuriteLengthDistribution)
            
            self.send_response(200)
            self.end_headers()
            self.wfile.write(b'{"status": "success"}')
        else:
            self.send_response(404)
            self.end_headers()
            self.wfile.write(b'{"status": "not found"}')

def run_server():
    with socketserver.TCPServer(("", PORT), APIHandler) as httpd:
        logger.info(f"Serving at port {PORT}")
        httpd.serve_forever()

if __name__ == "__main__":
    run_server()
