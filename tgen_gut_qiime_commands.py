import pandas as pd
import psycopg2
from pathlib import Path
from dotenv import load_dotenv
import os
import numpy as np
import re
from qiime2 import Artifact
import qiime2.plugins.demux.actions as demux_actions
import qiime2.plugins.dada2.actions as dada2_actions
from qiime2 import Metadata
import glob
import qiime2.plugins.metadata.actions as metadata_actions
import qiime2.plugins.feature_table.actions as feature_table_actions
import qiime2.plugins.feature_table.methods
#In the directory /labs/Microbiome/Tgen_CLIA_GMHI/python_qiime_api_analysis
  
def view_metadata():
  df = pd.read_csv('TGen_sample_metadata.tsv', sep='\t')
  print(df)
  sample_metadata_md = Metadata.load('TGen_sample_metadata.tsv')
  metadata_summ_viz, = metadata_actions.tabulate(
    input=sample_metadata_md,
  )
  metadata_summ_viz.save('metadata_summ.qzv')

def create_casava_demux_artifacts():
  #Get the file locations, get the newest files by looking at the miseq run numbers, create symbollic links in TGen_Gut_Casava_Fastqs, import the demuxed sequences into a qiime2 artifact
  env_path = '/labs/COVIDseq/COVIDpoint/covidpoint_updates/.env'
  load_dotenv(dotenv_path=env_path)
  token=os.environ["mesquite_pass"]
  conn = psycopg2.connect("host=mesquite.tgen.org dbname= gut_microbiome user=bvan-tassel password = "+token+" sslmode=disable gssencmode=disable")
  cursor = conn.cursor() 

  sql = """select fastq_r1_location, fastq_r2_location from sample_fastqs where fastq_format = 'casava'"""
  
  cursor.execute(sql)
  result = cursor.fetchall()
  #result = result[:2]
  r1_results, r2_results = map(list, zip(*result))
  results = r1_results + r2_results

  if not os.path.exists('TGen_Gut_Casava_Fastqs'):
    os.mkdir('TGen_Gut_Casava_Fastqs')
  
  for result in results:
    if not os.path.exists('TGen_Gut_Casava_Fastqs/' + result.split('/')[-2]):
      os.mkdir('TGen_Gut_Casava_Fastqs/' + result.split('/')[-2])
    if not os.path.exists('TGen_Gut_Casava_Fastqs/' + result.split('/')[-2] + '/' + result.split('/')[-1]):
      os.symlink(result, 'TGen_Gut_Casava_Fastqs/' + result.split('/')[-2] + '/' + result.split('/')[-1])

  runs = glob.glob('/labs/Microbiome/Tgen_CLIA_GMHI/python_qiime_api_analysis/TGen_Gut_Casava_Fastqs/*')

  if not os.path.exists('demux_artifacts_by_run'):
    os.mkdir('demux_artifacts_by_run')

  for run in runs:
    run = '/'.join(run.split('/')[-2:])
    demux_sequences = Artifact.import_data('SampleData[PairedEndSequencesWithQuality]', run, view_type='CasavaOneEightSingleLanePerSampleDirFmt')
    demux_sequences.save('demux_artifacts_by_run/' + run.split('/')[1] + '-casava-paired-end-demux.qza')

    demultiplexed_sequences_summ_viz, = demux_actions.summarize(
        data=demux_sequences,
    )
    demultiplexed_sequences_summ_viz.save('demux_artifacts_by_run/' + run.split('/')[1] + '-casava-paired-end-demux.qzv')

def denoise_casava_artifacts():
  
  #Denoise the separate runs and output the denoise artifacts to a folder.
  """
  sbatch --mem 10g qiime dada2 denoise-single \
  --p-trim-left 0 \
  --p-trunc-len 150 \
  --i-demultiplexed-seqs manifest-paired-end-demux.qza.qza \
  --o-representative-sequences manifest-rep-seqs.qza \
  --o-table manifest-table.qza \
  """
  #help(dada2)

  casava_demux_artifacts = os.listdir('demux_artifacts_by_run')
  #print(cas_demux_artifacts)

  if not os.path.exists('denoised_casava_artifacts'):
    os.mkdir('denoised_casava_artifacts')

  for demux_artifact in casava_demux_artifacts:
    #print(demux_artifact)
    prepend_name = demux_artifact.split('-casava')[0]

    demultiplexed_sequences_casava = Artifact.load('demux_artifacts_by_run/' + demux_artifact)

    cas_feature_table, cas_asv_sequences, dada2_stats = dada2_actions.denoise_paired(
      demultiplexed_seqs=demultiplexed_sequences_casava,
      trunc_len_f=150,
      trim_left_r=0,
      trunc_len_r=150,
      n_threads = 10,)

    if not os.path.exists('denoised_casava_artifacts/' + prepend_name):
      os.mkdir('denoised_casava_artifacts/' + prepend_name)

    cas_feature_table.save('denoised_casava_artifacts/' + prepend_name + '/' + prepend_name + '_cas_feature_table.qza')
    cas_asv_sequences.save('denoised_casava_artifacts/' + prepend_name + '/' + prepend_name +'_cas_asv_sequences.qza')
    dada2_stats.save('denoised_casava_artifacts/' + prepend_name + '/' + prepend_name + '_cas_dada2_stats.qza')

    stats_dada2_md_md = dada2_stats.view(Metadata)
    dada2_stats_summ_viz, = metadata_actions.tabulate(
      input=stats_dada2_md_md,
    )

    dada2_stats_summ_viz.save('denoised_casava_artifacts/' + prepend_name + '/' + prepend_name + '_cas_dada2_stats.qzv')

    sample_metadata_md = Metadata.load('TGen_sample_metadata.tsv')
    cas_feature_table_summ_viz, = feature_table_actions.summarize(
    table=cas_feature_table,
    sample_metadata=sample_metadata_md,
    )
    
    cas_asv_sequences_summ_viz, = feature_table_actions.tabulate_seqs(
        data=cas_asv_sequences,
    )
    cas_feature_table_summ_viz.save('denoised_casava_artifacts/' + prepend_name + '/' + prepend_name + '_cas_feature_table.qzv')
    cas_asv_sequences_summ_viz.save('denoised_casava_artifacts/' + prepend_name + '/' + prepend_name + '_cas_asv_sequences')

def merge_casava_denoise_artifacts():
  #help(feature_table)
  denoise_folders = glob.glob('/labs/Microbiome/Tgen_CLIA_GMHI/python_qiime_api_analysis/denoised_casava_artifacts/*')
  first_folder =  denoise_folders[0]
  first_table = glob.glob(first_folder + '/*feature*.qza')[0]

  merged_table = Artifact.load(first_table)
  denoise_folders = denoise_folders[1:]
  for folder in denoise_folders:
    table_to_merge_in = glob.glob(folder + '/*feature*')[0]
    table_to_merge_in = Artifact.load(table_to_merge_in)
    merged_table=feature_table_actions.merge(tables=[merged_table, table_to_merge_in])[0]
  #merged=feature_table.merge(tables=[table_1, table_2])
  #help(qiime2.plugins.feature_table.methods.merge)

  if not os.path.exists('final_merges'):
    os.mkdir('final_merges')

  merged_table.save('final_merges/casava_denoise_table_merged.qza')
  casava_denoise_table_merged_viz = feature_table_actions.summarize(table=merged_table)[0]
  casava_denoise_table_merged_viz.save('final_merges/casava_denoise_table_merged_viz.qzv')

def create_manifest_demux_artifacts():
  #Get the file locations, create a manifest file, manifest.csv, import the demuxed sequences into a qiime2 artifact
  env_path = '/labs/COVIDseq/COVIDpoint/covidpoint_updates/.env'
  load_dotenv(dotenv_path=env_path)
  token=os.environ["mesquite_pass"]
  conn = psycopg2.connect("host=mesquite.tgen.org dbname= gut_microbiome user=bvan-tassel password = "+token+" sslmode=disable gssencmode=disable")
  cursor = conn.cursor() 

  sql = "select fastq_r1_location, fastq_r2_location from sample_fastqs where fastq_format = 'manifest'"

  cursor.execute(sql)
  fastq_locations_results = cursor.fetchall()
  manifest_list = []
  i = 0
  for fastq_locations in fastq_locations_results:
    
    if len(fastq_locations) != 2:
      print('error')
      exit()
    r1_location = ''
    r2_location = ''
    for location in fastq_locations:
      if '_R1_' in location:
        r1_location = location
      elif '_R2_' in location:
        r2_location = location
      
    if r1_location == '' or r2_location == '':
      print('No r1 and r2')
      exit()
    manifest_list.append([i, r1_location, r2_location])
    i+=1
  manifest_df = pd.DataFrame(manifest_list, columns=['sample-id', 'forward-absolute-filepath', 'reverse-absolute-filepath'])

  if not os.path.exists('manifest_csvs'):
    os.mkdir('manifest_csvs')
  
  manifest_df['folder'] = manifest_df['forward-absolute-filepath'].apply(lambda x: x.split('/')[-2])
  manifest_set = manifest_df['folder'].unique()

  for manifest_folder_location in manifest_set:
    if not os.path.exists('manifest_csvs/' + manifest_folder_location):
      manifest_subset_df = manifest_df[manifest_df['folder'] == manifest_folder_location]
      manifest_subset_df.to_csv('manifest_csvs/' + manifest_folder_location + '.tsv', sep ='\t', index = False)
  
  manifest_files = [x + '.tsv' for x in manifest_set]

  if not os.path.exists('manifest_artifacts_by_run'):
    os.mkdir('manifest_artifacts_by_run')

  for manifest in manifest_files:
    print(manifest)
    demux_sequences = Artifact.import_data('SampleData[PairedEndSequencesWithQuality]', 'manifest_csvs/' + manifest, view_type='PairedEndFastqManifestPhred33V2')

    demux_sequences.save('manifest_artifacts_by_run/' + manifest.split('.')[0] + '_manifest-paired-end-demux.qza')
  

    demultiplexed_sequences_summ_viz, = demux_actions.summarize(
        data=demux_sequences,
    )
    demultiplexed_sequences_summ_viz.save('manifest_artifacts_by_run/' + manifest.split('.')[0] + '_manifest-paired-end-demux.qzv')

def denoise_manifest_artifacts():
  manifest_demux_artifacts = os.listdir('manifest_artifacts_by_run')
  if not os.path.exists('denoised_manifest_artifacts'):
    os.mkdir('denoised_manifest_artifacts')

  for demux_artifact in manifest_demux_artifacts:
    print(demux_artifact)
    prepend_name = demux_artifact.split('_manifest')[0]
    print(prepend_name)

    demultiplexed_sequences_manifest = Artifact.load('manifest_artifacts_by_run/' + demux_artifact)
    
    manifest_feature_table, manifest_asv_sequences, dada2_stats = dada2_actions.denoise_paired(
      demultiplexed_seqs=demultiplexed_sequences_manifest,
      trunc_len_f=150,
      trim_left_r=0,
      trunc_len_r=150,
      n_threads = 10,)

    if not os.path.exists('denoised_manifest_artifacts/' + prepend_name):
      os.mkdir('denoised_manifest_artifacts/' + prepend_name)

    manifest_feature_table.save('denoised_manifest_artifacts/' + prepend_name + '/' + prepend_name + '_manifest_feature_table.qza')
    manifest_asv_sequences.save('denoised_manifest_artifacts/' + prepend_name + '/' + prepend_name +'_manifest_asv_sequences.qza')
    dada2_stats.save('denoised_manifest_artifacts/' + prepend_name + '/' + prepend_name + '_manifest_dada2_stats.qza')

    stats_dada2_md_md = dada2_stats.view(Metadata)
    dada2_stats_summ_viz, = metadata_actions.tabulate(
      input=stats_dada2_md_md,
    )

    dada2_stats_summ_viz.save('denoised_manifest_artifacts/' + prepend_name + '/' + prepend_name + '_manifest_dada2_stats.qzv')

    sample_metadata_md = Metadata.load('TGen_Sample_metadata.tsv')
    manifest_feature_table_summ_viz, = feature_table_actions.summarize( table=manifest_feature_table, sample_metadata=sample_metadata_md,)
    manifest_asv_sequences_summ_viz, = feature_table_actions.tabulate_seqs(
        data=manifest_asv_sequences,
    )
    manifest_feature_table_summ_viz.save('denoised_manifest_artifacts/' + prepend_name + '/' + prepend_name + '_manifest_feature_table.qzv')
    manifest_asv_sequences_summ_viz.save('denoised_manifest_artifacts/' + prepend_name + '/' + prepend_name + '_manifest_asv_sequences.qzv')

def merge_manifest_casava_artifacts():
  casava_merged_table=Artifact.load('final_merges/casava_denoise_table_merged.qza')
  manifest_merged_table=Artifact.load('denoised_manifest_artifacts/PMI-TGN_Highlander_COH_CBM588_AlloHSCT_18193/PMI-TGN_Highlander_COH_CBM588_AlloHSCT_18193_manifest_feature_table.qza')
  merged_table=feature_table_actions.merge(tables=[casava_merged_table, manifest_merged_table])[0]
  merged_table.save('final_table_merged.qza')
  manifest_denoise_table_merged_viz = feature_table_actions.summarize(table=merged_table)[0]
  manifest_denoise_table_merged_viz.save('final_table_merged_viz.qzv')

if __name__ == "__main__":
  #view_metadata()

  #create_casava_demux_artifacts()
  #denoise_casava_artifacts()
  #merge_casava_denoise_artifacts()
  
  #create_manifest_demux_artifacts()
  #denoise_manifest_artifacts()

  merge_manifest_casava_artifacts()