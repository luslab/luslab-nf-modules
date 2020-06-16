    def run(file_in):
      file_out = file_in.replace('.bed', '_peaks.bed')
      df_in = pd.read_csv(file_in,
                          names = ["chrom", "start", "end", "name", "score", "strand"],
                          header=None, sep='\t')

      df_out = df_in[['chrom', 'strand', 'start', 'score']]

      df_out.sort_values(['chrom', 'strand', 'start'], ascending=[True, True, True], inplace=True)

      paraclu_input = file_in + '.paraclu_input'
      paraclu_output = file_in + '.paraclu_output'

      df_out.to_csv(paraclu_input, sep='\t', header=None, index=None)

      call(f'/home/aram/BI_tools/paraclu-9/paraclu 10 "{paraclu_input}" | /home/aram/BI_tools/paraclu-9/paraclu-cut.sh > "{paraclu_output}"', shell=True)
      df_in = pd.read_csv(paraclu_output,
                          names = ["sequence_name", "strand","start", "end", "number_of_positions",
                                  "sum_of_data_values", "min_density", "max_density"],
                          header=None, sep='\t')
      df_in['fourth_column'] = '.'
      df_out = df_in[['sequence_name', 'start', 'end', 'fourth_column', 'sum_of_data_values', 'strand']]
      df_out.sort_values(['sequence_name','start', 'end', 'strand'],
                        ascending=[True, True, True, True], inplace=True)
      df_out.to_csv(file_out, sep='\t', header=None, index=None)
      call(f'rm "{paraclu_input}"', shell=True)
      call(f'rm  "{paraclu_output}"', shell=True)