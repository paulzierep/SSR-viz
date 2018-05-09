# pssm_minimal
Minimal set-up for the difference PSSM command line tool

Example code:
./build_csv.py -i test_data/test_input.fasta -o test.csv -d -r '\|([a-zA-Z0-9\-_]+)$'

./plot_dpssm.py -i test_data/class_labels.csv -a test_data/temp_alignment.fasta -p 'test'