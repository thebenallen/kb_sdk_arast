#BEGIN_HEADER
import os
import sys
import shutil
import hashlib
import subprocess
import traceback
import uuid
import logging
import pprint
import json
import tempfile
import re
from datetime import datetime

import numpy as np

from Bio import SeqIO

from biokbase.workspace.client import Workspace as workspaceService


logging.basicConfig(format="[%(asctime)s %(levelname)s %(name)s] %(message)s",
                    level=logging.DEBUG)
logger = logging.getLogger(__name__)

#END_HEADER


class AssemblyRAST:
    '''
    Module Name:
    AssemblyRAST

    Module Description:
    A KBase module: AssemblyRAST

This sample module contains multiple assembly methods:

    run_kiki
    run_velvet
    run_miniasm
    run_spades
    run_idba
    run_megahit
    run_ray
    run_masurca
    run_a5
    run_a6

    '''

    ######## WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    #BEGIN_CLASS_HEADER
    workspaceURL = None

    # target is a list for collecting log messages
    def log(self, target, message):
        # we should do something better here...
        if target is not None:
            target.append(message)
        print(message)
        sys.stdout.flush()

    def create_temp_json(self, attrs):
        f = tempfile.NamedTemporaryFile(delete=False)
        outjson = f.name
        f.write(json.dumps(attrs))
        f.close()
        return outjson

    # combine multiple read library objects into a kbase_assembly_input
    def combine_read_libs(self, libs):
        pe_libs = []
        se_libs = []
        refs = []
        for libobj in libs:
            data = libobj['data']
            info = libobj['info']
            print(json.dumps(data))
            print(json.dumps(info))
            type_name = info[2].split('.')[1].split('-')[0]
            lib = dict()
            if type_name == 'PairedEndLibrary':
                if 'lib1' in data:
                    lib['handle_1'] = data['lib1']['file']
                elif 'handle_1' in data:
                    lib['handle_1'] = data['handle_1']
                if 'lib2' in data:
                    lib['handle_2'] = data['lib2']['file']
                elif 'handle_2' in data:
                    lib['handle_2'] = data['handle_2']
                if 'interleaved' in data:
                    lib['interleaved'] = data['interleaved']
                pe_libs.append(lib)
            elif type_name == 'SingleEndLibrary':
                if 'lib' in data:
                    lib['handle'] = data['lib']['file']
                elif 'handle' in data:
                    lib['handle'] = data['handle']
                se_libs.append(lib)

        assembly_input = { 'paired_end_libs': pe_libs,
                           'single_end_libs': se_libs,
                           'references': refs }
        logger.debug('kbase_assembly_input = {}'.format(json.dumps(assembly_input)))
        return assembly_input

    # template
    def arast_run(self, ctx, params, assembler='kiki'):
        output = None

        #### do some basic checks
        if 'workspace_name' not in params:
            raise ValueError('workspace_name parameter is required')
        if 'read_library_name' not in params:
            raise ValueError('read_library_name parameter is required')
        if 'output_contigset_name' not in params:
            raise ValueError('output_contigset_name parameter is required')
        min_contig_len = params.get('min_contig_len') or 300

        token = ctx['token']

        os.environ["KB_AUTH_TOKEN"] = token
        os.environ["ARAST_URL"] = '140.221.67.209' # testing on torino

        ws = workspaceService(self.workspaceURL, token=token)
        objects = ws.get_objects([{'ref': params['workspace_name']+'/'+params['read_library_name']}])

        libs = [objects[0]]
        wsid = objects[0]['info'][6]

        kbase_assembly_input = self.combine_read_libs(libs)
        tmp_data = self.create_temp_json(kbase_assembly_input)

        logger.info('Start {} assembler'.format(assembler))

        cmd = ['ar-run', '-a', assembler, '--data-json', tmp_data]
        logger.debug('CMD: {}'.format(' '.join(cmd)))

        p = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, shell=False)

        out, err = p.communicate()
        logger.debug(out)

        if p.returncode != 0:
            raise ValueError('Error running ar_run, return code: {}\n'.format(p.returncode))

        job_id = None
        match = re.search('(\d+)', out)
        if match:
            job_id = match.group(1)
        else:
            raise ValueError('No integer job ID found: {}\n'.format(out))

        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()*1000)
        output_dir = os.path.join(self.scratch, 'output.'+str(timestamp))
        output_contigs = os.path.join(output_dir, 'contigs.fa')
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        cmdstr = 'ar-get -j {} -w -p | ar-filter -l {} > {}'.format(job_id, min_contig_len, output_contigs)
        logger.debug('CMD: {}'.format(cmdstr))
        subprocess.check_call(cmdstr, shell=True)

        cmd = ['ar-get', '-j', job_id, '-w', '-r']
        logger.debug('CMD: {}'.format(' '.join(cmd)))
        ar_report = subprocess.check_output(cmd)

        cmd = ['ar-get', '-j', job_id, '-w', '-l']
        logger.debug('CMD: {}'.format(' '.join(cmd)))
        ar_log = subprocess.check_output(cmd)

        # Warning: this reads everything into memory!  Will not work if
        # the contigset is very large!
        contigset_data = {
            'id': '{}.contigset'.format(assembler),
            'source': 'User assembled contigs from reads in KBase',
            'source_id':'none',
            'md5': 'md5 of what? concat seq? concat md5s?',
            'contigs':[]
        }

        lengths = []
        for seq_record in SeqIO.parse(output_contigs, 'fasta'):
            contig = {
                'id': seq_record.id,
                'name': seq_record.name,
                'description': seq_record.description,
                'length': len(seq_record.seq),
                'sequence': str(seq_record.seq),
                'md5': hashlib.md5(str(seq_record.seq)).hexdigest()
            }
            lengths.append(contig['length'])
            contigset_data['contigs'].append(contig)


        provenance = [{}]
        if 'provenance' in ctx:
            provenance = ctx['provenance']
        # add additional info to provenance here, in this case the input data object reference
        provenance[0]['input_ws_objects']=[params['workspace_name']+'/'+params['read_library_name']]

        # save the contigset output
        new_obj_info = ws.save_objects({
                'id': wsid, # set the output workspace ID
                'objects':[
                    {
                        'type': 'KBaseGenomes.ContigSet',
                        'data': contigset_data,
                        'name': params['output_contigset_name'],
                        'meta': {},
                        'provenance': provenance
                    }
                ]
            })

        os.remove(tmp_data)
        shutil.rmtree(output_dir)

        # create a Report
        report = ''
        report += '============= Raw Contigs ============\n' + ar_report + '\n'

        report += '========== Filtered Contigs ==========\n'
        report += 'ContigSet saved to: '+params['workspace_name']+'/'+params['output_contigset_name']+'\n'
        report += 'Assembled into '+str(len(contigset_data['contigs'])) + ' contigs.\n'
        report += 'Average Length: '+str(sum(lengths)/float(len(lengths))) + ' bp.\n'

        # compute a simple contig length distribution
        bins = 10
        counts, edges = np.histogram(lengths, bins)
        report += 'Contig Length Distribution (# of contigs -- min to max basepairs):\n'
        for c in range(bins):
            report += '   '+str(counts[c]) + '\t--\t' + str(edges[c]) + ' to ' + str(edges[c+1]) + ' bp\n'

        print report

        reportObj = {
            'objects_created':[{'ref':params['workspace_name']+'/'+params['output_contigset_name'], 'description':'Assembled contigs'}],
            'text_message': report
        }

        reportName = '{}.report.{}'.format(assembler, job_id)
        report_obj_info = ws.save_objects({
                'id': wsid,
                'objects': [
                    {
                        'type': 'KBaseReport.Report',
                        'data': reportObj,
                        'name': reportName,
                        'meta': {},
                        'hidden': 1,
                        'provenance': provenance
                    }
                ]
            })[0]

        output = { 'report_name': reportName, 'report_ref': str(report_obj_info[6]) + '/' + str(report_obj_info[0]) + '/' + str(report_obj_info[4]) }

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method filter_contigs return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return output

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.workspaceURL = config['workspace-url']
        self.scratch = os.path.abspath(config['scratch'])
        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)
        #END_CONSTRUCTOR
        pass

    def run_kiki(self, ctx, params):
        # ctx is the context object
        # return variables are: output
        #BEGIN run_kiki
        output = self.arast_run(ctx, params, "kiki")
        #END run_kiki
        return [output]

    def run_velvet(self, ctx, params):
        # ctx is the context object
        # return variables are: output
        #BEGIN run_velvet
        output = self.arast_run(ctx, params, "velvet")
        #END run_velvet
        return [output]

    def run_miniasm(self, ctx, params):
        # ctx is the context object
        # return variables are: output
        #BEGIN run_miniasm
        output = self.arast_run(ctx, params, "miniasm")
        #END run_miniasm
        return [output]

    def run_spades(self, ctx, params):
        # ctx is the context object
        # return variables are: output
        #BEGIN run_spades
        output = self.arast_run(ctx, params, "spades")
        #END run_spades
        return [output]

    def run_idba(self, ctx, params):
        # ctx is the context object
        # return variables are: output
        #BEGIN run_idba
        output = self.arast_run(ctx, params, "idba")
        #END run_idba
        return [output]

    def run_megahit(self, ctx, params):
        # ctx is the context object
        # return variables are: output
        #BEGIN run_megahit
        output = self.arast_run(ctx, params, "megahit")
        #END run_megahit
        return [output]

    def run_ray(self, ctx, params):
        # ctx is the context object
        # return variables are: output
        #BEGIN run_ray
        output = self.arast_run(ctx, params, "ray")
        #END run_ray
        return [output]

    def run_masurca(self, ctx, params):
        # ctx is the context object
        # return variables are: output
        #BEGIN run_masurca
        output = self.arast_run(ctx, params, "masurca")
        #END run_masurca
        return [output]

    def run_a5(self, ctx, params):
        # ctx is the context object
        # return variables are: output
        #BEGIN run_a5
        output = self.arast_run(ctx, params, "a5")
        #END run_a5
        return [output]

    def run_a6(self, ctx, params):
        # ctx is the context object
        # return variables are: output
        #BEGIN run_a6
        output = self.arast_run(ctx, params, "a6")
        #END run_a6
        return [output]
