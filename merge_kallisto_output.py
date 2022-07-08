### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python merge_kallisto_output.py
					--in <INPUT_FOLDER>
					--gff <GFF_FILE>
					--tpms <TPM_OUTPUT_FILE>
					--counts <COUNTS_OUTPUT_FILE>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""


import os, sys, glob

# --- end of imports --- #

def load_counttable( counttable ):
	"""! @brief load data from counttable """
	
	counts = {}
	tpms = {}
	with open( counttable, "r" ) as f:
		f.readline()	#remove header
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			counts.update( { parts[0]: float( parts[3] ) } )
			tpms.update( { parts[0]: float( parts[4] ) } )
			line = f.readline()
	return counts, tpms


def generate_mapping_table( gff_file ):
	"""! @brief generate transcript to gene mapping table """
	
	transcript2gene = {}
	with open( gff_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] in [ "mRNA", "transcript" ]:
					try:
						ID = parts[-1].split(';')[0].split('=')[1]
						parent = parts[-1].split('arent=')[1]
						if ";" in parent:
							parent = parent.split(';')[0]
						transcript2gene.update( { ID: parent } )
					except:
						print line
			line = f.readline()
	return transcript2gene


def map_counts_to_genes( transcript2gene, counts ):
	"""! @brief map transcript counts to parent genes """
	
	error_collector = []
	gene_counts = {}
	for key in counts.keys():
		try:
			gene_counts[ transcript2gene[ key ] ] += counts[ key ]
		except KeyError:
			try:
				gene_counts.update( { transcript2gene[ key ]: counts[ key ] } )
			except KeyError:
				error_collector.append( key )
				#print key
				gene_counts.update( { key: counts[ key ] } )
	if len( error_collector ) > 0:
		print "number of unmapped transcripts: " + str( len( error_collector ) )
	return gene_counts


def generate_output_file( output_file, data ):
	"""! @brief generate output file for given data dictionary """
	
	samples = sorted( data.keys() )
	
	with open( output_file, "w" ) as out:
		out.write( "\t".join( [ "gene" ] + samples ) + '\n' )
		for gene in sorted( data.values()[0].keys() ):
			new_line = [ gene ]
			for sample in samples:
				new_line.append( data[ sample ][ gene ] )
			out.write( "\t".join( map( str, new_line ) ) + '\n' )


def main( arguments ):
	"""! @brief run everything """
	
	data_input_dir = arguments[ arguments.index( '--in' )+1 ]
	gff_file = arguments[ arguments.index( '--gff' )+1 ]
	if '--counts' in arguments:
		counts_output_file = arguments[ arguments.index( '--counts' )+1 ]
	else:
		counts_output_file = False
	if '--tpms' in arguments:
		tpm_output_file = arguments[ arguments.index( '--tpms' )+1 ]
	else:
		tpm_output_file = False
	
	transcript2gene = generate_mapping_table( gff_file )
	print "number of mapped transcripts: " + str( len( transcript2gene.keys() ) )
	counttables = glob.glob( data_input_dir + "*.tsv" )
	print "number of detected counttables: " + str( len( counttables ) )

	count_data = {}
	tpm_data = {}
	for filename in counttables:
		ID = filename.split('/')[-1].split('.')[0]
		counts, tpms = load_counttable( filename )
		#TPM are available and could be processed in the same way
		gene_counts = map_counts_to_genes( transcript2gene, counts )
		gene_tpms = map_counts_to_genes( transcript2gene, tpms )
		count_data.update( { ID: gene_counts } )
		tpm_data.update( { ID: gene_tpms } )
	
	if counts_output_file:
		generate_output_file( counts_output_file, count_data )
	if tpm_output_file:
		generate_output_file( tpm_output_file, tpm_data )


if '--in' in sys.argv and '--gff' in sys.argv and '--tpms' in sys.argv:
	main( sys.argv )
elif '--in' in sys.argv and '--gff' in sys.argv and '--counts' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
