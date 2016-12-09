#!/usr/bin/python
from SimpleHTTPServer import SimpleHTTPRequestHandler
import BaseHTTPServer
import pdb

import dnds

class CORSRequestHandler (SimpleHTTPRequestHandler):
    # def end_headers (self):
    #     self.send_header('Access-Control-Allow-Origin', '*')
    #     self.send_header('Access-Control-Allow-Methods', 'GET, POST, PUT, DELETE, OPTIONS')
    #     self.send_header('Access-Control-Allow-Headers', '*')
    #     SimpleHTTPRequestHandler.end_headers(self)

    def _set_headers(self):
        self.send_response(200)
        self.send_header('Content-type', 'text/html')
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Methods', 'GET, POST, PUT, DELETE, OPTIONS')
        self.send_header('Access-Control-Allow-Headers', '*')
        self.end_headers()

    def do_OPTIONS(self):           
        self.send_response(200, "ok")       
        self.send_header('Access-Control-Allow-Origin', '*')                
        self.send_header('Access-Control-Allow-Methods', 'GET, POST, OPTIONS')
        self.send_header("Access-Control-Allow-Headers", "X-Requested-With")        

    def do_GET(self):           
        self.send_response(200)
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Content-type',    'text/html')                                    
        self.end_headers()              
        self.wfile.write("<html><body>Hello world!</body></html>")
        self.connection.shutdown(1) 

    def do_POST(self):
        content_length = int(self.headers['Content-Length']) # <--- Gets the size of data
        post_data = self.rfile.read(content_length) # <--- Gets the data itself

        qry_seq_raw, ref_seq_raw = post_data.split("+")

        dnds_data = dnds.dnds_pipeline(qry_seq_raw, ref_seq_raw)

        #pdb.set_trace()

        self._set_headers()
        self.wfile.write("<html><body><h1>"+str(dnds_data)+"</h1></body></html>")

if __name__ == '__main__':
    BaseHTTPServer.test(CORSRequestHandler, BaseHTTPServer.HTTPServer)