import logging
import json
import os
import pprint
import requests
import subprocess
import sys
import time
import unittest

from os import environ
from ConfigParser import ConfigParser
from requests_toolbelt import MultipartEncoder
from pprint import pprint

from biokbase.workspace.client import Workspace as workspaceService
from biokbase.AbstractHandle.Client import AbstractHandle as HandleService
from AssemblyRAST.AssemblyRASTImpl import AssemblyRAST

logging.basicConfig(format="[%(asctime)s %(levelname)s %(name)s] %(message)s",
                    level=logging.INFO)
logger = logging.getLogger(__name__)

class AssemblyRASTTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        cls.ctx = {'token': token, 'provenance': [{'service': 'AssemblyRAST',
            'method': 'please_never_use_it_in_production', 'method_params': []}],
            'authenticated': 1}
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('AssemblyRAST'):
            cls.cfg[nameval[0]] = nameval[1]
        cls.wsURL = cls.cfg['workspace-url']
        cls.ws = workspaceService(cls.wsURL, token=token)
        cls.serviceImpl = AssemblyRAST(cls.cfg)

        cls.shockURL = cls.cfg['shock-url']
        cls.handleURL = cls.cfg['handle-service-url']
        logger.info('shockURL='+cls.shockURL)


    @classmethod
    def tearDownClass(cls):
        logger.info('tearDownClass() temporarily disabled')
        # if hasattr(cls, 'wsName'):
        #     cls.ws.delete_workspace({'workspace': cls.wsName})
        #     print('Test workspace was deleted')


    def getWsClient(self):
        return self.__class__.ws


    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_AssemblyRAST_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})
        self.__class__.wsName = wsName
        return wsName


    def getImpl(self):
        return self.__class__.serviceImpl


    def getContext(self):
        return self.__class__.ctx


    # call this method to get the WS object info of a Paired End Library (will
    # upload the example data if this is the first time the method is called during tests)
    def getPairedEndLibInfo(self):

        if hasattr(self.__class__, 'pairedEndLibInfo'):
            return self.__class__.pairedEndLibInfo

        # try to reuse persist test json file if it exists
        testFile = 'data/testPairedEndLibInfo.json'
        if os.path.exists(testFile):
            logger.info("Reading pairedEndLibInfo from {}".format(testFile))
            with open(testFile) as testInfoFile:
                return json.load(testInfoFile)

        # 1) upload files to shock
        token = self.ctx['token']
        forward_shock_file = self.curl_upload_file_to_shock(
            shock_service_url = self.shockURL,
            filePath = 'data/small.forward.fq',
            token = token
            )
        reverse_shock_file = self.curl_upload_file_to_shock(
            shock_service_url = self.shockURL,
            filePath = 'data/small.reverse.fq',
            token = token
            )
        #pprint(forward_shock_file)
        #pprint(reverse_shock_file)

        # 2) create handle
        hs = HandleService(url=self.handleURL, token=token)
        forward_handle = hs.persist_handle({
                                        'id' : forward_shock_file['id'],
                                        'type' : 'shock',
                                        'url' : self.shockURL,
                                        'file_name': forward_shock_file['file']['name'],
                                        'remote_md5': forward_shock_file['file']['checksum']['md5']})

        reverse_handle = hs.persist_handle({
                                        'id' : reverse_shock_file['id'],
                                        'type' : 'shock',
                                        'url' : self.shockURL,
                                        'file_name': reverse_shock_file['file']['name'],
                                        'remote_md5': reverse_shock_file['file']['checksum']['md5']})

        # 3) save to WS
        paired_end_library = {
            'lib1': {
                'file': {
                    'hid':forward_handle,
                    'file_name': forward_shock_file['file']['name'],
                    'id': forward_shock_file['id'],
                    'url': self.shockURL,
                    'type':'shock',
                    'remote_md5':forward_shock_file['file']['checksum']['md5']
                },
                'encoding':'UTF8',
                'type':'fastq',
                'size':forward_shock_file['file']['size']
            },
            'lib2': {
                'file': {
                    'hid':reverse_handle,
                    'file_name': reverse_shock_file['file']['name'],
                    'id': reverse_shock_file['id'],
                    'url': self.shockURL,
                    'type':'shock',
                    'remote_md5':reverse_shock_file['file']['checksum']['md5']
                },
                'encoding':'UTF8',
                'type':'fastq',
                'size':reverse_shock_file['file']['size']

            },
            'interleaved':0,
            'sequencing_tech':'artificial reads'
        }

        new_obj_info = self.ws.save_objects({
                        'workspace':self.getWsName(),
                        'objects':[
                            {
                                'type':'KBaseFile.PairedEndLibrary',
                                'data':paired_end_library,
                                'name':'test.pe.reads',
                                'meta':{},
                                'provenance':[
                                    {
                                        'service':'AssemblyRAST',
                                        'method':'test_kiki'
                                    }
                                ]
                            }]
                        })
        self.__class__.pairedEndLibInfo = new_obj_info[0]

        logger.info("pairedEndLibInfo='{}'".format(json.dumps(new_obj_info[0])))

        return new_obj_info[0]


    # this method somehow fails with a 'bad connection' error, use the curl_ version instead
    # Helper script borrowed from the transform service, logger removed
    def upload_file_to_shock(self,
                             shock_service_url = None,
                             filePath = None,
                             ssl_verify = True,
                             token = None):
        """
        Use HTTP multi-part POST to save a file to a SHOCK instance.
        """

        if token is None:
            raise Exception("Authentication token required!")

        #build the header
        header = dict()
        header["Authorization"] = "Oauth {0}".format(token)

        if filePath is None:
            raise Exception("No file given for upload to SHOCK!")

        dataFile = open(os.path.abspath(filePath), 'rb')
        m = MultipartEncoder(fields={'upload': (os.path.split(filePath)[-1], dataFile)})
        header['Content-Type'] = m.content_type

        logger.info("Sending {0} to {1}".format(filePath,shock_service_url))
        try:
            response = requests.post(shock_service_url + "/node", headers=header, data=m, allow_redirects=True, verify=ssl_verify)
            dataFile.close()
        except:
            dataFile.close()
            raise

        if not response.ok:
            response.raise_for_status()

        result = response.json()

        if result['error']:
            raise Exception(result['error'][0])
        else:
            return result["data"]

    # borrowed from the assembly service
    def curl_upload_file_to_shock(self,
                             shock_service_url = None,
                             filePath = None,
                             ssl_verify = True,
                             token = None):

        cmd = ['curl',
               '-X', 'POST',
               '-F', 'upload=@{}'.format(filePath),
               '{}/node/'.format(shock_service_url)]

        cmd += ['-s']
        cmd += ['-H', '"Authorization: OAuth {}"'.format(token)]

        sys.stderr.write("Uploading: {}\n".format(' '.join(cmd)))
        logger.debug("curl_post_file: {}".format(' '.join(cmd)))
        r = subprocess.check_output(' '.join(cmd), shell=True)
        sys.stderr.write("\n")

        res = json.loads(r)
        logger.debug("res: {}\n".format(res))
        # pprint(res)

        if res['error']:
            raise Exception(res['error'][0])
        else:
            return res["data"]



    # Example test method for filtering contigs
    # def test_filter_contigs_ok(self):
    #     obj_name = "contigset.1"
    #     contig1 = {'id': '1', 'length': 10, 'md5': 'md5', 'sequence': 'agcttttcat'}
    #     contig2 = {'id': '2', 'length': 5, 'md5': 'md5', 'sequence': 'agctt'}
    #     contig3 = {'id': '3', 'length': 12, 'md5': 'md5', 'sequence': 'agcttttcatgg'}
    #     obj1 = {'contigs': [contig1, contig2, contig3], 'id': 'id', 'md5': 'md5', 'name': 'name',
    #             'source': 'source', 'source_id': 'source_id', 'type': 'type'}
    #     self.getWsClient().save_objects({'workspace': self.getWsName(), 'objects':
    #         [{'type': 'KBaseGenomes.ContigSet', 'name': obj_name, 'data': obj1}]})
    #     ret = self.getImpl().filter_contigs(self.getContext(), {'workspace': self.getWsName(),
    #         'contigset_id': obj_name, 'min_length': '10'})
    #     obj2 = self.getWsClient().get_objects([{'ref': self.getWsName()+'/'+obj_name}])[0]['data']
    #     self.assertEqual(len(obj2['contigs']), 2)
    #     self.assertTrue(len(obj2['contigs'][0]['sequence']) >= 10)
    #     self.assertTrue(len(obj2['contigs'][1]['sequence']) >= 10)
    #     self.assertEqual(ret[0]['n_initial_contigs'], 3)
    #     self.assertEqual(ret[0]['n_contigs_removed'], 1)
    #     self.assertEqual(ret[0]['n_contigs_remaining'], 2)


    def test_run_kiki(self):

        # figure out where the test data lives
        pe_lib_info = self.getPairedEndLibInfo()
        pprint(pe_lib_info)

        # Object Info Contents
        # 0 - obj_id objid
        # 1 - obj_name name
        # 2 - type_string type
        # 3 - timestamp save_date
        # 4 - int version
        # 5 - username saved_by
        # 6 - ws_id wsid
        # 7 - ws_name workspace
        # 8 - string chsum
        # 9 - int size
        # 10 - usermeta meta

        params = {
            'workspace_name': pe_lib_info[7],
            'read_library_name': pe_lib_info[1],
            'output_contigset_name': 'output.contigset',
            'min_contig_length': 350,
            #'extra_params': '-k 23'
        }

        result = self.getImpl().run_kiki(self.getContext(),params)
        print('RESULT:')
        pprint(result)
