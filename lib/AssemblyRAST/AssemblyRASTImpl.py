#BEGIN_HEADER
import os
import sys
import traceback
import uuid

from biokbase.workspace.client import Workspace as workspaceService
#END_HEADER


class AssemblyRAST:
    '''
    Module Name:
    AssemblyRAST

    Module Description:
    A KBase module: AssemblyRAST
This sample module contains one small method - filter_contigs.
    '''

    ######## WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    #BEGIN_CLASS_HEADER
    workspaceURL = None
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.workspaceURL = config['workspace-url']
        #END_CONSTRUCTOR
        pass

    def run_kiki(self, ctx, params):
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN run_kiki

        print('Starting run_kiki method.')

        # if 'workspace' not in params:
        #     raise ValueError('Parameter workspace is not set in input arguments')
        # workspace_name = params['workspace']
        # try:
        #     min_length = int(min_length_orig)
        # except ValueError:
        #     raise ValueError('Cannot parse integer from min_length parameter (' + str(min_length_orig) + ')')
        # if min_length < 0:
        #     raise ValueError('min_length parameter shouldn\'t be negative (' + str(min_length) + ')')

        token = ctx['token']
        ws = workspaceService(self.workspaceURL, token=token)

        objects = ws.get_objects([{'ref': params['workspace_name']+'/'+params['read_library_name']}])

        data = objects[0]['data']
        info = objects[0]['info']
        type_name = info[2].split('.')[1].split('-')[0]

        report = 'report will go here\n'
        report += '\tinput data type: '+type_name
        reportObj = {
            # 'objects_created':[{'ref':params['workspace_name']+'/'+params['output_contigset_name'], 'description':'Assembled contigs'}],
            'objects_created':[],
            'text_message':report
        }


        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects']=[params['workspace_name']+'/'+params['read_library_name']]

        reportName = 'megahit_report_'+str(hex(uuid.getnode()))
        report_obj_info = ws.save_objects({
                'id':info[6],
                'objects':[
                    {
                        'type':'KBaseReport.Report',
                        'data':reportObj,
                        'name':reportName,
                        'meta':{},
                        'hidden':1,
                        'provenance':provenance
                    }
                ]
            })[0]

        output = { 'report_name': reportName, 'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]) }

        #END run_kiki

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method filter_contigs return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [output]
